import os
import sys
import time
import glob
import shutil
import subprocess
from lib.tools.job_launcher import launch_jobs
from lib.parsing.dir_parsing_tools import generate_fastq_table


def generate_annotation_fasta(gtf, genome, outdir, programs):

    if not genome:
        print(time.asctime(), f"ERROR: Genome FASTA is required to generate transcriptome FASTA file")
        return None

    print(time.asctime(), f"Generating transcriptome FASTA of annotation: {gtf}")

    gtf_outname = os.path.basename(gtf).replace(".gtf", "")

    t_fasta = os.path.join(outdir, f"{gtf_outname}.fa")
    gffread_exe = programs["gffread"]

    command = f"{gffread_exe} -w {t_fasta} -g {genome} {gtf}"

    print(time.asctime(), "GffRead command:")
    print(command + "\n")

    subprocess.call(command, shell=True)

    # print(time.asctime(), f"Transcripts fasta file: {t_fasta}")

    return t_fasta


def generate_salmon_index(trans_fasta, outdir, programs):

    print(time.asctime(), f"Generating Salmon index file")

    t_index = os.path.join(outdir, "transcripts_index")

    if os.path.exists(t_index):
        print(time.asctime(), f'WARNING: Salmon transcriptome index "{t_index}" already exist')
        print(time.asctime(), f'Using current index\n')
        return t_index

    # Only available option for Salmon +v.1.2.0 (previous version used fmd option)
    index_type = "quasi"  # Recommended in older versions
    index_type = "puff"   # Only option in latest version

    salmon_exe = programs["salmon"]
    command = f"{salmon_exe} index -t {trans_fasta} -i {t_index} --type {index_type} --keepDuplicates"

    print(time.asctime(), f"Salmon index creation command:")
    print(command + "\n")

    subprocess.call(command, shell=True)

    # print(time.asctime(), f"Transcripts index file: {t_index}")

    return t_index


def generate_salmon_commands(fastq_table, trans_index, outdir, programs, paired=True, n_threads=2):

    print("\n")
    print(time.asctime(), f"Generating Salmon quantification commands")

    salmon_exe = programs["salmon"]

    if not paired:
        sys.exit(f'ERROR: The current version if only able to handle paired-end data.')

    commands_list = []
    # Get paired files from table
    with open(fastq_table) as fh:
        # Skip header
        next(fh)
        for row in fh:
            fl_1, fl_2 = row.strip("\n").split(",")

            if not fl_2:
                sys.exit(f'ERROR: The current version if only able to handle paired-end data.')

            # Use the name of file 1 as the outname for the analysis
            outname = os.path.splitext(os.path.basename(fl_1))[0]
            outfile = os.path.join(outdir, outname)

            # This tell Salmon to automatically detect the library type
            lib_type = "A"

            command = f"{salmon_exe} quant -i {trans_index} -l {lib_type} -1 {fl_1} -2 {fl_2} -o {outfile} --useVBOpt --seqBias --threads {n_threads}"

            commands_list.append(command)

    return commands_list


def launch_salmon_analysis(gtf, genome_fl, fastq_table, paths_dt, paired=True):

    # Get absolute path
    salmon_exe = shutil.which("salmon")
    gffread_exe = shutil.which("gffread")

    if not salmon_exe:
        print(time.asctime(), "ERROR: Salmon is not installed in your environment.")

    if not gffread_exe:
        print(time.asctime(), "ERROR: Gffread is not installed in your environment.")

    programs = {"salmon": salmon_exe, "gffread": gffread_exe}

    # Generate transcriptome fasta file
    trans_fa = generate_annotation_fasta(gtf, genome_fl, paths_dt["quant"], programs)

    # Generate Salmon index
    trans_index = generate_salmon_index(trans_fa, paths_dt["quant"], programs)

    # Generate Salmon quantification command for each sample
    salmon_commands = generate_salmon_commands(fastq_table, trans_index, paths_dt["quant"], programs, paired=paired, n_threads=4)

    test_log = os.path.join(paths_dt["logs"], "salmon_commands.txt")
    with open(test_log, "w+") as fh:
        for command in salmon_commands:
            fh.write(command + "\n\n")

    print("\n")
    print(time.asctime(), f"Logfile of Salmon quantification commands: {test_log}\n")

    # Use multiprocessing to run as many parallel quantification as possible
    launch_jobs(salmon_commands, log_dir=paths_dt["logs"])

    return paths_dt["quant"]


def quantify_transcripts(gtf, genome_fl, fastq_dir, paths_dt, tool="Salmon"):

    print(time.asctime(), f"Starting transcripts quantification")

    # Write FASTQ files into a table, useful to track how the files were paired for the analysis
    fastq_table = generate_fastq_table(fastq_dir, paths_dt["report"])

    # Quantify transcripts
    if tool.upper() == "SALMON":
        quant_dir = launch_salmon_analysis(gtf, genome_fl, fastq_table, paths_dt, paired=True)
    else:
        sys.exit(f"ERROR: Quantification tool {tool} not supported")

    # Check that the number of output quantification output files match the number of input files
    ext = "sf"  # Extension of Salmon output
    quant_files = glob.glob(os.path.join(quant_dir, f"**/*.{ext}"), recursive=True)
    if not quant_files:
        sys.exit(f'ERROR: The directory containing Salmon quantification output is empty.')

    fastq_files = []
    with open(fastq_table) as fh:
        next(fh)  # Skip header
        for row in fh:
            fastq_files.append(row)

    if not len(quant_files) == len(fastq_files):
        sys.exit(f"ERROR: There are Salmon output files missing in the quantification folder: {quant_dir}")

    return quant_dir

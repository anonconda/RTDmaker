import os
import sys
import time
import glob
import warnings
from collections import defaultdict

from lib.report.report_tools import simple_write_table, percentiles_from_counts
from lib.parsing.gtf_object_tools import create_gtf_object

warnings.filterwarnings("ignore")


def generate_sj_report(sj_dir, outpath, outname):

    print(time.asctime(), f"Generating report of SJ dataset", flush=True)

    # Create outfolder if it doesn't exit
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    outfile = os.path.join(outpath, outname)

    # 1) Set up some values for the analysis

    # Require number of uniquely mapped reads
    sj_read_th = 1

    # Ignore reads without enough overhang, make this argument optional?
    overhang_th = 10

    # 2) Get the SJ dataset files

    # Get SJ files generated by STAR
    sj_files = sorted(glob.glob(os.path.join(sj_dir, "**/*SJ.out.tab"), recursive=True))

    # Remove 1st pass log files
    sj_files = [fl for fl in sj_files if "STARpass1" not in fl]

    # Check if folder is empty
    if not sj_files:
        sys.exit(f'ERROR: Folder {sj_dir} does not contain any "SJ.out.tab" files.')

    # 3) Create dictionary to track for each SJ co-ordinates the N of reads and the file

    # Format: SJ_ID: (N reads, Filename). Ex: 'Chr01-:100-200': (10, 'File01_SJ.out.tab')
    sj_reads_dt = defaultdict(list)

    # 4) At the same time, create dictionary with read counts to generate a percentile distribution table per each file
    file_count_dt = {}

    # Track the assigned filenames to fill missing values in the dictionary further downstream
    observed_files = set()

    for i, sj_file in enumerate(sorted(sj_files, key=os.path.getmtime), 1):

        sj_name = os.path.basename(sj_file)
        sj_name = f"{i:03}_{sj_name}"

        observed_files.add(sj_name)

        # This dictionary is to track the counts of uniquely mapped reads in the file. This dict is then stored
        sj_count_dt = defaultdict(int)

        with open(sj_file) as fh:
            for row in fh:
                row_lst = row.strip("\n").split("\t")
                # SJ.out.tab - high confidence collapsed splice junctions in tab-delimited format.
                # Column 0: chromosome
                # Column 1: first base of the intron (1-based)
                # Column 2: last base of the intron (1-based)
                # Column 3: strand (0:  undefined, 1:  +, 2:  -)
                # Column 4: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
                # Column 5: 0: unannotated, 1: annotated (only if splice junctions database is used)
                # Column 6: number of uniquely mapping reads crossing the junction
                # Column 7: number of multi-mapping reads crossing the junction
                # Column 8: maximum spliced alignment overhang

                if row_lst[3] == '1':
                    strand = '+'
                elif row_lst[3] == '2':
                    strand = '-'
                else:
                    strand = '.'

                # Adding strand information to make sure it is unique
                sj_id = f"{row_lst[0]}{strand}:{row_lst[1]}-{row_lst[2]}"

                # Report only SJ with enough overhang
                check_overhang = True
                if check_overhang:
                    if int(row_lst[8]) < overhang_th:
                        continue

                # Check if SJ is canonical
                check_canonical = True
                if check_canonical:
                    intron_motif = int(row_lst[4])
                    # Canonical "flags" according to STAR documentation
                    if intron_motif not in (1, 2, 3, 4, 5, 6):
                        continue

                # Number of uniquely mapped reads to SJ
                n_unique_reads = int(row_lst[6])

                # Store information for the SJ-reads table
                sj_reads_dt[sj_id].append((n_unique_reads, sj_name))

                # Increase read count for the per-file percentile table
                sj_count_dt[n_unique_reads] += 1

        # Store dict for per-file percentile table
        file_count_dt[sj_name] = sj_count_dt

    # Fill missing values in SJ-reads dict (i.e. SJ not present in a file should report value 0 for that file name
    updated_sj_reads_dt = defaultdict(list)
    for sj_id, sj_data in sj_reads_dt.items():

        updated_data = sj_data

        sj_reads, sj_files = list(zip(*sj_data))

        for file_tag in observed_files.difference(sj_files):
            updated_data.append((0, file_tag))

        # IMPORTANT: elements MUST have the same order across all values to write the table, this order is used for head
        updated_sj_reads_dt[sj_id] = sorted(updated_data, key=lambda tpl: tpl[1])

    header = "SJ_ID\t" + "\t".join(sorted(observed_files)) + "\n"
    sj_rows = [header]

    for sj_id, sj_data in updated_sj_reads_dt.items():

        sj_reads, sj_files = list(zip(*sj_data))

        row = f"{sj_id}\t" + "\t".join([str(v) for v in sj_reads]) + "\n"
        sj_rows.append(row)

    # Write SJ-ID - Reads table
    simple_write_table(sj_rows, f"{outfile}_SJ_reads_numbers.tsv")

    # Write SJ per-file percentile distribution table
    perc = [0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99, 100]

    header = "File\t" + "\t".join([str(v) for v in perc]) + "\n"
    sj_perc_rows = [header]
    for sj_file, count_dt in file_count_dt.items():
        count_dt = percentiles_from_counts(count_dt, percentiles_range=perc)
        sj_quants, sj_vals = list(zip(*count_dt))

        row = f"{sj_file}\t" + "\t".join([str(v) for v in sj_vals]) + "\n"
        sj_perc_rows.append(row)

    simple_write_table(sj_perc_rows, f"{outfile}_SJ_reads_percentiles.tsv")

    print(time.asctime(), f"Report tables location: {outpath}", flush=True)

    return outpath

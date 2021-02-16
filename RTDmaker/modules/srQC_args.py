"""
Created on Mon Jul 31 07:07:07 2017
@author: Juan Carlos Entizne
@email: e.entizne[at]dundee.ac.uk
"""

import os
import sys
import time
import traceback
import modules.srQC_main as srQC
from lib.tools.logger import logger, clean_log
from lib.tools.input_tools import check_input, create_project_dir
from argparse import ArgumentParser, RawTextHelpFormatter


description = \
    "Description:\n" + \
    "RTDmaker 'ShortReads' module process transcriptome annotations assembled from RNA-seq data to generate a " \
    "High-Quality RTD annotation.\n" \
    "This module identifies and removes:\n" \
    "1) Redundant Transcripts, 2) Transcripts with vague location information, " \
    "3) Transcripts poorly supported by reads (at SJ and overall expression level), " \
    "4) Fragmentary Transcripts, and 5) Chimeric transcripts.\n" \
    "The output of this module is a single transriptome annotation, RTD, " \
    "containing only high-quality and well supported transcripts models."

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)

parser.add_argument('--assemblies',
                    dest="ass_dir",
                    action="store", default=None,
                    help="Path of the folder containing the transcriptome assemblies (GTF format) to be analyzed.")

parser.add_argument('--references',
                    dest="ref_dir",
                    action="store", default=None,
                    help="Path of the folder containing reference annotations (GTF format) to be integrated into the RTD.")

###

parser.add_argument('--SJ-data',
                    dest="sj_dir",
                    action="store", default=None,
                    help="Path of the folder containing the splice-junction data ('SJ.out.tab' files) generated by STAR.")

parser.add_argument("--SJ-reads",
                    dest="sjreads_th", default=(2, 1),
                    nargs=2, metavar=("n_reads", "n_samples"),
                    help="Minimum number of uniquely-map reads (X) observed in a minimum number of samples (Y) to consider "
                         "a splice-junction (SJ) supported. Default: 2 1 (2 uniquely-map reads in at least 1 sample).")

###

parser.add_argument('--genome',
                    dest="genome_fl",
                    action="store", default=None,
                    help="Path of the Genome FASTA file to extract the transcripts sequence to perform transcript quantification.")

parser.add_argument('--fastq',
                    dest="fastq_dir",
                    action="store", default=None,
                    help="Path of the folder containing the FASTQ files to perform transcript quantification. Supported extensions: 'fq.gz', 'fq', 'fastq.gz', 'fastq'.")

parser.add_argument("--tpm",
                    dest="abund_th",
                    default=(1, 1),
                    nargs=2, metavar=("n_tpm", "n_samples"),
                    help="Minimum TPM abundance (N), observed in a minimum number of samples (M), to consider a transcript "
                         "as supported by quantification. Default: 1 1 (1 TPM in at least 1 sample).")

###

parser.add_argument("--fragment-len",
                    dest="len_th",
                    type=float, default=0.7,
                    help="Minimum percentage of gene-length coverage below which a transcript is identified as "
                         "fragment. Value must be within 0 and 1. Default: 0.7")

parser.add_argument("--antisense-len",
                    dest="antlen_th",
                    type=float, default=0.5,
                    help="Minimum percentage of gene-length coverage below which a 'monoexonic-antisense' transcript "
                         "in the opposite strand is identified as fragment. Value must be within 0 and 1. Default: 0.5")

###

parser.add_argument("--add",
                    dest="add_set", nargs='*',
                    default=['None'], choices=['unstranded', 'intronic', 'None'],
                    help="Potentially non-coding transcripts to include into the annotation. "
                         "Options: 'unstranded', 'intronic'. Default: exclude them.")

###

parser.add_argument("--ram",
                    dest="size_th",
                    type=int, default=8,
                    help="Maximum size (in Gb) of temporary files. Change this value only if your system cannot handle "
                         "uploading into memory files with the current file size. Default: 8 Gb.")

parser.add_argument("--keep",
                    dest="keep_set", nargs='*',
                    default=['None'], choices=['intermediary', 'removed', 'None'],
                    help="Intermediary files to keep during the analysis. "
                         "Options: 'intermediary', 'removed'. Default: Remove them.")

###

parser.add_argument('--outpath',
                    dest="outpath", default=None,
                    help="Path of the output folder.")

parser.add_argument('--prefix',
                    dest="prefix", default=None,
                    help="Prefix to assign to the genes and transcripts IDs in the RTD. Default: outname.")

parser.add_argument('--outname',
                    dest="outname", default=None,
                    help="Name of the output file.")


def main():

    args = parser.parse_args()

    # The pipeline allowed initially to optionally avoid doing the SJ and/or Abundance QC via the arg 'args.skip_set'
    # However, these QC analysis are too important, so I have decided to disabled the original option here
    # To avoid breaking the code, I just pass the original argument empty
    args.skip_set = set()

    # Check user input
    args = check_input(args)

    # Create project folder and subfolder structure
    paths_dt = create_project_dir(args.outpath)

    # Create logfile to track the analysis
    time_stamp = time.strftime("%Y%m%d-%H%M%S")

    logfile = os.path.join(paths_dt["logs"], f"{time_stamp}_{args.outname}_logfile_temporary.txt")
    logger(logfile, w_mode="w+")

    print("\n")
    print(time.asctime(), "Starting RTDmaker Short-Reads analysis\n")

    # Record command typed by the user
    command = " ".join(sys.argv)
    print("Command:")
    print(f"{command}\n", flush=True)

    # Unhandled exception - bad news. We want to record the traceback to our log, not just stderr
    try:

        srQC.main(args.ass_dir, args.ref_dir,                       # Annotations to analyze
                  args.sj_dir, args.sjreads_th,                     # SJ-QC related args
                  args.genome_fl, args.fastq_dir, args.abund_th,    # Abundance-QC related args
                  args.len_th, args.antlen_th,                      # Structure-QC related args
                  args.skip_set, args.add_set,                      # Category-selection related args
                  args.size_th, args.keep_set,                      # Memory-related arguments
                  paths_dt, args.prefix, args.outname, logfile)     # Output related args

    except SystemExit as err:
        # Valid for python 3.5+
        print("".join(traceback.TracebackException.from_exception(err).format()))
    except Exception as err:
        print("".join(traceback.TracebackException.from_exception(err).format()))
        sys.exit(f"{err}")

    # Remove redundant lines from logfile
    clean_log(logfile)
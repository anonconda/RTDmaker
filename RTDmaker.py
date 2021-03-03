"""
Created on Mon Jul 31 07:07:07 2017
@author: Juan Carlos Entizne
@email: e.entizne[at]dundee.ac.uk
"""

import sys
import argparse
from logging import getLogger
from modules import srQC_args as srQC  # short-reads Quality-Control (srQC)

description = \
    "Description:\n" + \
    "RTDmaker is a computational pipeline to generate High-Quality transcriptome annotations, known as Reference-Transcript-Datasets (RTDs). " \
    "Currently, RTDmaker is made of one module, 'ShortReads', that process transcriptome annotations assembled from RNA-seq data. " \
    "Contact: Juan.Carlos.Entizne@hutton.ac.uk\n"

parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter, prog='RTDmaker')

# Version
parser.add_argument('-v', '--version', action='version', version='0.1.5')
subparsers = parser.add_subparsers()

# Parser of the module for RNA-seq analysis
ShortReads_Subparser = subparsers.add_parser("ShortReads", parents=[srQC.parser],
                                             help="Quality Control of transcriptome assemblies from RNA-seq data.")
ShortReads_Subparser.set_defaults(which="ShortReads")

# Setting logging preferences
logger = getLogger(__name__)


def main():

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    try:
        args = parser.parse_args()
        # Setting modules parser
        if args.which == "ShortReads":
            srQC.parser = parser
            srQC.main()
        else:
            print(f"Command {args.which} not recognized.")
            sys.exit(1)

    except Exception as e:
        logger.error(f"{e}")
        sys.exit(1)


if __name__ == "__main__":
    main()


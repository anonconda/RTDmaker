import os
import sys
import time
import glob
import argparse

parser = argparse.ArgumentParser(description="Create a table contained the FASTQ files in a directory")

parser.add_argument('--fastq-dir', dest="fastq_dir", action="store", default=None,
                    help="Path of directory containing the FASTQ files")

parser.add_argument("--tags", dest="tags_list", nargs='*', default=None,
                    help="Text that must present in the file name to keep it in the table. Default: None (keep all)")

parser.add_argument('--paired', dest="paired", default=False, action='store_true',
                    help="Boolean. Whether the FASTQ files in the directory are paired or not")

parser.add_argument('--outpath', dest="outpath", default=None, help="Path of the output file")


def generate_fastq_table(fastq_dir, outpath, file_tags=None, paired=True):

    # print(time.asctime(), f"Generating table reporting FASTQ files in directory: {fastq_dir}")

    if not file_tags:
        file_tags = set()

    for ext in ['fq.gz', 'fq', 'fastq.gz', 'fastq']:
        temp, observed = [], set()
        for filename in sorted(glob.glob(os.path.join(fastq_dir, f"**/*.{ext}"), recursive=True)):

            if file_tags:
                for tag in file_tags:
                    if tag in filename:
                        temp.append(filename)
                        observed.add(ext)
            else:
                temp.append(filename)
                observed.add(ext)
        temp.sort()

        # This will stop the execution once a type of extension is add into the table
        if observed and temp:
            break

    if not temp:
        sys.exit(f"No FastQ files found in the specified directory: {fastq_dir}")

    if not paired:
        sys.exit(f"ERROR: Currently this script support only the paired FASTQ data")

    paired_files = []
    if paired:
        for fl_pair in [temp[i:i + 2] for i in range(0, len(temp), 2)]:
            fl_pair = tuple(fl_pair)
            paired_files.append(fl_pair)
    else:
        for fl in temp:
            fl_pair = (fl, "")
            paired_files.append(fl_pair)

    table = os.path.join(outpath, "fastq_table.csv")
    with open(table, "w+") as fh:
        if paired:
            fh.write("File_pair_1,File_pair_2\n")
        else:
            fh.write("File\n")
        for file_pair in paired_files:
            fh.write(f"{file_pair[0]},{file_pair[1]}\n")

    print(time.asctime(), f"FASTQ table location: {table}")

    return table


if __name__ == "__main__":
    args = parser.parse_args()
    generate_fastq_table(args.fastq_dir, args.outpath, file_tags=args.tags_list, paired=args.paired)

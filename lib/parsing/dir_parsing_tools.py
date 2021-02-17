import os
import sys
import time
import glob


def generate_annotation_table(dir_path, outpath, outname):

    if not dir_path:
        return None

    # Currently, the RTDmaker pipeline only process annotations in GTF format
    ext = "gtf"
    temp = []
    for filename in sorted(glob.glob(os.path.join(dir_path, f"**/*.{ext}"), recursive=True)):

        # Ignore Cufflinks ignored models annotations
        if filename.endswith("skipped.gtf"):
            continue
        else:
            temp.append(filename)
    temp.sort()

    if not temp:
        sys.exit(f"No GTF files found in the specified directory: {dir_path}")

    # Write table
    table = os.path.join(outpath, f"{outname}_table.csv")
    with open(table, "w+") as fh:
        fh.write("File_path,File_tag\n")
        for i, fl_path in enumerate(temp, 1):
            fh.write(f"{fl_path},file{i}\n")

    return table


def parse_annotation_table(input_table, force=False):

    # print(time.asctime(), f'Parsing input table: {input_table}')

    try:
        gtf_files, tags = [], set()
        with open(input_table) as fh:
            next(fh)
            for row in fh:
                file_path, file_tag = row.strip("\n").split(",")
                file_path = file_path.replace('"', '')
                gtf_files.append((file_path, file_tag))
                tags.add(file_tag)
    except Exception as err:
        print(err)
        sys.exit(f'ERROR while parsing input table {input_table}')

    if len(tags) != len(gtf_files):
        if force:
            print("WARNING: Number of tags does not match the number of annotation files in the input table")
        else:
            sys.exit(f"ERROR: Number of tags does not match number of annotation files. Please check input table.")

    return gtf_files


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

import os
import sys
import time
import glob
from collections import defaultdict
from lib.report.report_tools import global_track_dt
from lib.parsing.gtf_object_tools import create_gtf_object, write_gtf


def remove_low_abundant_via_quantification(gtf_file, quant_dir, paths_dt, outname, quant_th=None, to_keep=None):

    if not to_keep:
        to_keep = set()

    print(time.asctime(), "Starting abundance Quality Control analysis")

    # Check and sanitize threshold input
    if not quant_th:
        quant_th = (0, 0)
    assert len(quant_th) == 2

    ext = "sf"  # Extension of Salmon output
    quant_files = glob.glob(os.path.join(quant_dir, f"**/*.{ext}"), recursive=True)
    if not quant_files:
        sys.exit(f'ERROR: The directory containing Salmon quantification output is empty.')

    # Unpack the threshold values
    try:
        tpm_th, sample_th = [float(val) for val in quant_th]
    except ValueError as err:
        sys.exit(err)

    ext = "sf"  # File extension of Salmon quantification output
    quant_dt = defaultdict(list)
    for filename in sorted(glob.glob(os.path.join(quant_dir, f"**/*.{ext}"), recursive=True)):

        # print(time.asctime(), f"Parsing quantification file: {filename}")

        with open(filename) as fh:
            # Skip header: Name	Length	EffectiveLength	TPM	NumReads
            next(fh)
            for row in fh:
                t_id, t_len, t_efflen, t_tpm, t_reads = row.strip("\n").split("\t")

                t_tpm = float(t_tpm)
                quant_dt[t_id].append(t_tpm)

    print(time.asctime(), "Identifying low-abundant transcripts")
    abundant, low_abundant = [set() for _ in range(2)]
    for t_id, t_abundances in sorted(quant_dt.items()):

        # Identify high-abundant values (TPM above selected TPM threshold abundance)
        high_abundance = [t_tpm for t_tpm in t_abundances if t_tpm >= tpm_th]

        # Identify if the abundance is present in enough sampels
        if len(high_abundance) >= sample_th:
            abundant.add(t_id)
        else:
            low_abundant.add(t_id)

    # print(time.asctime(), f"Removing low-abundant transcripts")
    gtf_obj = create_gtf_object(gtf_file)

    outfile = write_gtf(gtf_obj, abundant, paths_dt["inter"], f"{outname}_abundant.gtf")

    if 'removed' in to_keep:
        _ = write_gtf(gtf_obj, low_abundant, paths_dt["removed"], f"{outname}_low_abundant.gtf")

    # Track numbers of accepted / removed
    global_track_dt[f"Removed#{outname}_low_abundant.gtf"] = len(low_abundant)
    global_track_dt[f"Accepted#{outname}_abundant.gtf"] = len(abundant)

    return outfile

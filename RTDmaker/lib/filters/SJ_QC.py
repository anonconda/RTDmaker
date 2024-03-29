import os
import sys
import time
import glob
import warnings
import linecache
from collections import defaultdict
from lib.report.report_tools import global_track_dt
from lib.parsing.gtf_object_tools import create_gtf_object

warnings.filterwarnings("ignore")


def identify_supported_sj(sj_dir, sjreads_th, paths_dt, outname, verb=False):

    print(time.asctime(), "Identifying supported Splice Junctions", flush=True)

    # Read overhang value
    overhang_th = 10

    # Unpack the threshold values
    try:
        reads_th, samples_th = [int(val) for val in sjreads_th]
    except ValueError as err:
        sys.exit(err)

    # Get SJ files generated by STAR
    sj_files = sorted(glob.glob(os.path.join(sj_dir, "**/*SJ.out.tab"), recursive=True))

    # Remove 1st pass log files
    # print(time.asctime(), f'Ignoring "SJ.out.tab" files from 1st STAR pass (contain "STARpass1" in filename)')
    sj_files = [fl for fl in sj_files if "STARpass1" not in fl]

    # Check if folder is empty
    if not sj_files:
        sys.exit(f'ERROR: Folder {sj_dir} does not contain any "SJ.out.tab" files.')

    # Dict tracking the number of samples that show enough support for a SJ
    sj_support_dt = defaultdict(int)

    # Classification of rejected SJ
    non_canonical, low_reads = [set() for _ in range(2)]

    # print(time.asctime(), f'Identifying supported Splice Junctions from STAR "SJ.out.tab" files')
    for i, sj_file in enumerate(sj_files, 1):

        # sj_name = os.path.basename(sj_file)
        # print(time.asctime(), f'Parsing STAR Splice Junction file: {sj_file} ({i} of {len(sj_files)})')

        with open(sj_file) as fh:
            for row in fh:
                row_l = row.strip("\n").split("\t")
                # SJ.out.tab - high confidence collapsed splice junctions in tab-delimited format.
                # Only junctions supported by uniquely mapping reads are reported.
                # Column 0: chromosome
                # Column 1: first base of the intron (1-based)
                # Column 2: last base of the intron (1-based)
                # Column 3: strand (0:  undefined, 1:  +, 2:  -)
                # Column 4: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
                # Column 5: 0: unannotated, 1: annotated (only if splice junctions database is used)
                # Column 6: number of uniquely mapping reads crossing the junction
                # Column 7: number of multi-mapping reads crossing the junction
                # Column 8: maximum spliced alignment overhang

                # A SJ is accepted if:
                # a) It is non-canonical
                # b) It show enough reads (*) in enough examples
                # c) (*) the reads should have enough overhand

                # Adding strand information to make sure it is unique
                if row_l[3] == '1':
                    strand = '+'
                elif row_l[3] == '2':
                    strand = '-'
                else:
                    strand = '.'
                sj_id = f"{row_l[0]}{strand}:{row_l[1]}-{row_l[2]}"

                # c) Not enough overhand, ignore
                if int(row_l[8]) < overhang_th:
                    continue

                # b) Not canonical, ignore
                intron_motif = int(row_l[4])
                # Canonical "flags" according to STAR documentation
                if intron_motif in (1, 2, 3, 4, 5, 6):
                    pass
                else:
                    non_canonical.add(sj_id)
                    # By using continue here we are completely disregarding the read support of non-canonical SJ
                    continue

                # a) If it have enough reads, increase the number of samples that support the SJ
                n_unique_reads = int(row_l[6])
                if n_unique_reads >= reads_th:
                    sj_support_dt[sj_id] += 1

    # If the memory consumption is too high, I can explore changing this dictionary into a cPickle files instead
    # print(time.asctime(), f'Selecting supported Splice Junctions', flush=True)

    # Dictionary containing ONLY supported SJ
    supported_sj_dt = {}
    for sj_id, n_supported in sj_support_dt.items():
        if n_supported >= samples_th:
            supported_sj_dt[sj_id] = n_supported
        else:
            if sj_id not in non_canonical:
                low_reads.add(sj_id)

    # Write a table containing the rejected SJ
    # print(time.asctime(), f"Generating rejected SJ tables")
    for rejected_set, tag in zip([non_canonical, low_reads], ["non_canonical", "insufficient_reads"]):
        # Avoid writing tables if it is empty
        if rejected_set:
            rejected_outfile = os.path.join(paths_dt['removed'], f"{outname}_table_{tag}_SJ.csv")

            # print(time.asctime(), f'Writing output file: {rejected_outfile}')
            with open(rejected_outfile, "w+") as fh:
                fh.write("SJ_id\n")
                for sj_id in sorted(rejected_set):
                    # Remove strand information from the sj_id to make it compatible with IGV visualization tool
                    split_id = sj_id.split(":")
                    sj_id = f"{split_id[0][:-1]}:{split_id[1]}"
                    fh.write(f"{sj_id}\n")

    return supported_sj_dt


def remove_unsupported_SJ_models(gtf_file, supported_sj_dt, paths_dt, outname, to_keep=None):

    if not to_keep:
        to_keep = set()

    print(time.asctime(), f"Removing transcripts containing poorly supported Splice Junctions", flush=True)

    gtf_obj = create_gtf_object(gtf_file)

    print(time.asctime(), f'Identifying transcripts with supported Splice Junctions')
    supported_transcripts, monoexons = [set() for _ in range(2)]
    for trans_id, trans_introns in gtf_obj.trans_introns_dt.items():
        chrom_strand = gtf_obj.trans_chrom_dt[trans_id]

        # Important! Accept trans_id without introns (monoexons) as they can be valid but they can't show support here
        is_supported = True
        if not trans_introns:
            monoexons.add(trans_id)
            continue

        for intron in trans_introns:
            # Recreate SJ id
            sj_key = f"{chrom_strand}:{intron[0]}-{intron[1]}"

            # If one SJ fails its support, remove the whole transcript away
            try:
                _ = supported_sj_dt[sj_key]
            except KeyError:
                # If it is not in the dictionary (and it is not a monoexon), then it is not supported, remove it
                is_supported = False

        if is_supported:
            supported_transcripts.add(trans_id)

    # Re-introduce mono-exonic transcripts into the accepted output
    supported_transcripts = supported_transcripts | monoexons
    non_supported_trans = set(gtf_obj.trans_exons_dt.keys()) - supported_transcripts

    supported_out = os.path.join(paths_dt['inter'], f"{outname}_supported_SJ.gtf")
    nonsupported_out = os.path.join(paths_dt['removed'], f"{outname}_unsupported_SJ.gtf")

    # print(time.asctime(), f'Writing output file: {supported_out}')
    # print(time.asctime(), f'Writing output file: {nonsupported_out}')

    # Not using the write_file() method in gtf_object_tools as it is more efficient (memory and speed wise) in this case
    # to use the if/else approach to write into the proper file
    with open(supported_out, "w+") as fa, open(nonsupported_out, "w+") as fb:
        for trans_id in gtf_obj.trans_gene_dt.keys():
            for line_ix in gtf_obj.trans_gtf_lines_index[trans_id]:
                line = linecache.getline(gtf_obj.gtf_path, line_ix)
                if trans_id in supported_transcripts:
                    fa.write(line)
                else:
                    if 'removed' in to_keep:
                        fb.write(line)

    # Remove empty file
    if 'removed' not in to_keep:
        os.remove(nonsupported_out)

    # Track numbers of accepted / removed
    global_track_dt[f"Removed#{outname}_unsupported_SJ.gtf"] = len(non_supported_trans)
    global_track_dt[f"Accepted#{outname}_supported_SJ.gtf"] = len(supported_transcripts)

    return supported_out

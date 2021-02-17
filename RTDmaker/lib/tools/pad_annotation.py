import os
import sys
import time
import linecache
from collections import defaultdict

from lib.tools.logger import logger
from lib.parsing.gtf_object_tools import create_gtf_object
from lib.tools.other_tools import group_transcripts_by_overlap, sort_by_start, flat


def pad_transcripts(strand_trans, gtf_obj):

    strand_trans = sort_by_start(gtf_obj, strand_trans)
    overlap_groups = group_transcripts_by_overlap(gtf_obj, strand_trans, strict=False)

    modified_exons_dt = {}
    for subgroup in overlap_groups:
        # Find smallest and largest exon-coordinates in this subgroup
        subgroup_exons = []
        for trans_id in subgroup:
            subgroup_exons.extend(flat(gtf_obj.trans_exons_dt[trans_id]))
        min_exon, max_exon = min(subgroup_exons), max(subgroup_exons)

        # Substitute the first and last exon of the transcripts with the min/max values
        for trans_id in subgroup:
            prev_trans_exons = gtf_obj.trans_exons_dt[trans_id]

            if len(prev_trans_exons) < 2:
                modified_trans_exons = [(min_exon, max_exon)]

            else:
                mid_exons = prev_trans_exons[1:-1]
                first_exon = (min_exon, prev_trans_exons[0][-1])
                last_exon = (prev_trans_exons[-1][0], max_exon)
                modified_trans_exons = [first_exon] + mid_exons + [last_exon]

            modified_exons_dt[trans_id] = modified_trans_exons

    assert len(strand_trans) == len(modified_exons_dt.keys())

    return modified_exons_dt


def write_modified_gtf(gtf_obj, transcripts, feature, outfolder, outname, w_mode="w+", all=False):

    # For now handle only "exon" features
    if feature != "exon":
        sys.exit('Error. For the moment this method handle only features of type "exon"')

    if feature == "exon":
        feature_dt = gtf_obj.trans_exons_dt
    else:
        sys.exit('Error. For the moment this method handle only features of type "exon"')

    # Sort transcripts by the key (Chrom, Gene_ID, Trans_ID, Leftmost_coordinate)
    get_sort_key = lambda t_id: (gtf_obj.trans_chrom_dt[t_id], gtf_obj.trans_gene_dt[t_id], t_id,
                                 gtf_obj.trans_exons_dt[t_id][0][0])

    if all:
        # print(time.asctime(), "Writing all transcripts in input annotation into the output file")
        transcripts = set(gtf_obj.trans_exons_dt.keys())

    sorted_transcripts = sorted(transcripts, key=lambda t_id: get_sort_key(t_id))

    outfile = os.path.join(outfolder, outname)
    # print(time.asctime(), f'Writing output file: {outfile}')
    with open(outfile, w_mode) as fh:
        for trans_id in sorted_transcripts:
            if trans_id in transcripts:
                # Get representative line for the transcript
                trans_ix = gtf_obj.trans_gtf_lines_index[trans_id][0]
                rep_line = linecache.getline(gtf_obj.gtf_path, trans_ix)
                try:
                    seqname, source, _, _, _, score, strand, frame, attr = rep_line.split("\t")
                except ValueError:
                    sys.exit(f"Error. Line {trans_ix} of file {gtf_obj.gtf_path} seems to be broken:\n{rep_line}\n")
                score = "."

                # Select the new features
                for (start, end) in feature_dt[trans_id]:
                    mod_line = f"{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attr}\n"
                    fh.write(mod_line)

                for line_ix in gtf_obj.trans_gtf_lines_index[trans_id]:
                    line = linecache.getline(gtf_obj.gtf_path, line_ix)
                    # If the line feature (3rd columnd of GTF files) is the one modified, avoid re-annotation
                    line_feature = line.split("\t")[2]
                    if line_feature != feature:
                        fh.write(line)

    return outfile


def pad_annotation(gtf_file, paths_dt, outname, logfile):

    logger(logfile)

    print(time.asctime(), f'Generating padded annotation for file: {gtf_file}', flush=True)

    gtf_obj = create_gtf_object(gtf_file)

    # Group Trans by their strand of origin
    chrom_strand_transcripts_dt = defaultdict(list)
    for trans_id, _ in gtf_obj.trans_exons_dt.items():
        chrom_strand = gtf_obj.trans_chrom_dt[trans_id]
        chrom_strand_transcripts_dt[chrom_strand].append(trans_id)

    # print(time.asctime(), "Identifying gene boundaries")
    # First, group trans by their overlap.
    # Then, find smallest/largest exon-coord in the group and re-annotated it to all transcripts in overlap group.
    global_modified_exons_dt = {}
    for chrom_strand, strand_trans in sorted(chrom_strand_transcripts_dt.items()):
        modified_exons_dt = pad_transcripts(strand_trans, gtf_obj)
        global_modified_exons_dt.update(modified_exons_dt)

    print(time.asctime(), "Padding transcripts coordinates")
    # Re-annotate the transcripts exons in the gtf_obj
    gtf_obj = gtf_obj._replace(trans_exons_dt=global_modified_exons_dt)

    # I can't use the method write_gtf() as I have to modify the exon coordinates
    padded_gtf = write_modified_gtf(gtf_obj, set(), "exon", paths_dt["outpath"], f"{outname}_padded.gtf", all=True)

    print(time.asctime(), f"Padding completed: {padded_gtf}")

    return padded_gtf

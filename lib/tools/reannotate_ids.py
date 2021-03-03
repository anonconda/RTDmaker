import os
import time
import linecache
from collections import defaultdict

from lib.tools.logger import logger
from lib.report.report_tools import global_track_dt
from lib.parsing.gtf_object_tools import create_gtf_object
from lib.tools.other_tools import group_transcripts_by_overlap, flat


def get_group_start(trans_group, gtf_obj):

    group_exons = set()
    for t_id in trans_group:
        t_exons = gtf_obj.trans_exons_dt[t_id]
        group_exons.update(flat(t_exons))

    return min(group_exons)


def reannotate_ids(gtf_file, prefix, outfolder, outname):

    print(time.asctime(), f"Re-annotating genes and transcript IDs", flush=True)

    gtf_obj = create_gtf_object(gtf_file)

    # print(time.asctime(), "Generating novel IDs")
    # Dictionaries to track Old_id - New_id to generate a lookup table downstream
    gene_reannotation_dt, transcript_reannotation_dt = [{} for _ in range(2)]

    # Create a re-annotation gene dictionary
    i, step = 1, 10
    processed_genes, processed_trans = [set() for _ in range(2)]

    # 1st) We group transcripts by their overlap at the strand level
    # Note: Cause we want the Gene ID numbering to be done across strand we save the overlapping groups into a new dict
    chrom_overlap_groups_dt = defaultdict(list)
    for chrom, strand_transcripts in sorted(gtf_obj.chrom_trans_dt.items()):

        # Use un-stranded chromosome/scaffold as key (We want the Gene numbering to be incremental across strands)
        if not chrom[-1] in {"+", "-", "."}:
            print(time.asctime(), f'WARNING: Scaffold "{chrom}" does not finish with recognized strand symbol.')

        chrom_key = chrom[:-1]

        # Group transcripts by their overlap, use refine=True to avoid grouping transcripts with just small overlaps
        overlapping_trans = group_transcripts_by_overlap(gtf_obj, strand_transcripts, feature="exon", strict=False)

        chrom_overlap_groups_dt[chrom_key].extend(overlapping_trans)

    # Assign a Gene ID number to each overlapping group
    for chrom_key, overlapping_trans in chrom_overlap_groups_dt.items():

        # Sort overlap groups according to their leftmost exon (to establish Gene numbering from smallest to largest)
        overlapping_trans = sorted(overlapping_trans, key=lambda t_group: get_group_start(t_group, gtf_obj))

        for i, overlap_group in enumerate(overlapping_trans, 1):
            # Create novel Gene ID for each overlapping group
            novel_gene_id = f"{prefix}_{chrom_key}G{i * step:06}"

            trans_n = 1
            for prev_trans_id in sorted(overlap_group, key=lambda t_id: gtf_obj.trans_exons_dt[t_id][0]):

                prev_gene_id = gtf_obj.trans_gene_dt[prev_trans_id]
                if prev_gene_id not in processed_genes:
                    gene_reannotation_dt[prev_gene_id] = novel_gene_id
                    processed_genes.add(prev_gene_id)

                # Create novel Transcript ID
                if prev_trans_id not in processed_trans:
                    novel_trans_id = f"{novel_gene_id}.RTD.{trans_n}"
                    transcript_reannotation_dt[prev_trans_id] = novel_trans_id
                    processed_trans.add(prev_trans_id)
                    trans_n += 1

    # print(time.asctime(), "Re-annotating to novel IDs")
    final_rtd = f"{outname}.gtf"
    outfile = os.path.join(outfolder, final_rtd)
    with open(outfile, "w+") as fh:
        for chrom, gene_list in gtf_obj.chrom_gene_dt.items():
            for prev_gene_id in gene_list:
                for prev_trans_id in gtf_obj.gene_trans_dt[prev_gene_id]:
                    for line_ix in gtf_obj.trans_gtf_lines_index[prev_trans_id]:
                        prev_line = linecache.getline(gtf_obj.gtf_path, line_ix)

                        novel_gene_id = gene_reannotation_dt[prev_gene_id]
                        novel_trans_id = transcript_reannotation_dt[prev_trans_id]

                        # Remove the extra information from the attribute field
                        pre_attr_field = prev_line.split("\t")[-1]
                        new_attr_field = f'transcript_id "{novel_trans_id}"; gene_id "{novel_gene_id}";\n'

                        # temp_1 = prev_line.replace(pre_attr_field, new_attr_field)
                        # temp_2 = temp_1.replace(prev_gene_id, novel_gene_id)

                        novel_line = prev_line.replace(pre_attr_field, new_attr_field)

                        fh.write(novel_line)

    # Track numbers of accepted / removed
    global_track_dt[f"Final#{outname}.gtf"] = len(gtf_obj.trans_exons_dt.keys())

    return outfile, (gene_reannotation_dt, transcript_reannotation_dt)


def write_id_lookup_table(reannotation_dicts, outfolder, outname):

    # print(time.asctime(), f"Generating Gene/Transcript IDs look-up table")

    # Create output subfolder if it doesn't exist
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    gene_reannotation_dt, transcript_reannotation_dt = reannotation_dicts

    # Write table tracking the re-annotations
    for set_tag, reannotation_dt in zip(["Genes", "Transcripts"], [gene_reannotation_dt, transcript_reannotation_dt]):
        table_path = os.path.join(outfolder, f"{outname}_table_ID_lookup_{set_tag}.csv")
        with open(table_path, "w+") as fh:
            # Header
            fh.write("Original_ID,Novel_id\n")
            for (old_id, novel_id) in sorted(reannotation_dt.items()):
                fh.write(f"{old_id},{novel_id}\n")


def generate_reannotated_gtf(gtf_file, prefix, paths_dt, outname, logfile):

    logger(logfile)

    reannotated_gtf, reannotation_dicts = reannotate_ids(gtf_file, prefix, paths_dt["outpath"], outname)

    lookup_subfolder = os.path.join(paths_dt["report"], "ID_lookup_tables")
    write_id_lookup_table(reannotation_dicts, lookup_subfolder, outname)

    # print(time.asctime(), f"Gene/Transcript IDs re-annotation completed: {reannotated_gtf}")

    return reannotated_gtf


def generate_mapping_table(gtf_file, paths_dt, outname):

    print(time.asctime(), f"Generating genes to transcripts mapping table", flush=True)

    gtf_obj = create_gtf_object(gtf_file)

    outfile = os.path.join(paths_dt["report"], f"{outname}_table_ID_mapping.csv")
    with open(outfile, "w+") as fh:
        fh.write("TXNAME,GENEID\n")
        for gene_id, trans_list in sorted(gtf_obj.gene_trans_dt.items()):
            for trans_id in sorted(trans_list):
                fh.write(f"{trans_id},{gene_id}\n")

    print(time.asctime(), f"Mapping table completed: {outfile}", flush=True)

    return outfile

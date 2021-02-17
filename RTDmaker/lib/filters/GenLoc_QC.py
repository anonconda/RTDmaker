import os
import time
from collections import defaultdict
from itertools import combinations, permutations

from lib.report.report_tools import global_track_dt
from lib.filters.structural_QC import identify_intronic_transcripts
from lib.parsing.gtf_object_tools import create_gtf_object, write_gtf
from lib.tools.other_tools import get_transcript_length, get_longest_transcript
from lib.tools.other_tools import group_transcripts_by_overlap, group_transcripts_across_strands


def identify_unstranded_transcripts(gtf_obj):

    print(time.asctime(), f"Identifying unstranded transcripts")

    invalid_strand_trans = set()
    for trans_id, trans_strand in gtf_obj.trans_sense_dt.items():
        if trans_strand not in {"+", "-"}:
            invalid_strand_trans.add(trans_id)

    return invalid_strand_trans


def identify_monoexon_gene_transcripts(gtf_obj, trans_group):

    overlapping_transcripts = group_transcripts_by_overlap(gtf_obj, trans_group)
    single_monoexonic = set()
    for overlap_group in overlapping_transcripts:
        # Only transcript in group
        if len(overlap_group) == 1:
            trans_id = overlap_group[0]
            trans_exons = gtf_obj.trans_exons_dt[trans_id]
            # If it is mono-exonic
            if len(trans_exons) == 1:
                single_monoexonic.add(trans_id)

    return single_monoexonic


def identify_monoexonic_antisense(gtf_obj, to_ignore=None):

    if to_ignore is None:
        to_ignore = set()

    print(time.asctime(), f"Identifying mono-exonic antisense transcripts")

    # print(time.asctime(), f"Grouping chromosomes/scaffolds")
    # Group related strands; the strands can be "+", "-", and "." (unknown)
    grouped_chrom_dt = defaultdict(list)
    for chrom, trans_list in gtf_obj.chrom_trans_dt.items():
        if not chrom[-1] in {"+", "-", "."}:
            print(f'WARNING: Chromosome/Scaffold "{chrom}" does not finish with strand tag. Please check.')

        # Remove strand tag
        chrom_key = chrom[:-1]
        grouped_chrom_dt[chrom_key].append(chrom)

    # An antisense monoexon can be either completly nested withing the gene in the other strand, or just partially overlap it
    strict_nested = False

    # Relate gene coordinates to transcript IDs to track against which models is the selected transcripts an antisense
    gene_coord_trans_dt, antisense_relation_dt = [defaultdict(set) for _ in range(2)]

    # print(time.asctime(), f"Analyzing overlap of Mono-exonic Genes across strands")
    monoexonic_antisense_set = set()
    for chrom_key, chrom_groups in grouped_chrom_dt.items():
        # combinations() is faster but it doesn't allow to reveal all the possible relationship of antisense transcripts
        for chrom_A, chrom_B in permutations(chrom_groups, 2):

            strand_A_trans = gtf_obj.chrom_trans_dt[chrom_A]
            strand_B_trans = gtf_obj.chrom_trans_dt[chrom_B]

            strand_A_trans = set(t_id for t_id in strand_A_trans if t_id not in to_ignore)
            strand_B_trans = set(t_id for t_id in strand_B_trans if t_id not in to_ignore)

            if not strand_A_trans or not strand_B_trans:
                continue

            # print(time.asctime(), f'Processing scaffolds: {chrom_A} / {chrom_B}')

            # Get information from the 1st chromosome
            strand_A_monoexonics = identify_monoexon_gene_transcripts(gtf_obj, strand_A_trans)

            strand_A_gene_coordinates = set()
            for gene_id in gtf_obj.chrom_gene_dt[chrom_A]:
                gene_c = gtf_obj.gene_coords_dt[gene_id]
                strand_A_gene_coordinates.add(gene_c)

                gene_transcripts = gtf_obj.gene_trans_dt[gene_id]
                gene_coord_trans_dt[tuple(gene_c)].update(gene_transcripts)

            strand_A_gene_coordinates = sorted(strand_A_gene_coordinates)

            # Get information from the 2nd chromosome
            strand_B_monoexonics = identify_monoexon_gene_transcripts(gtf_obj, strand_B_trans)

            strand_B_gene_coordinates = set()
            for gene_id in gtf_obj.chrom_gene_dt[chrom_B]:
                gene_c = gtf_obj.gene_coords_dt[gene_id]
                strand_B_gene_coordinates.add(gene_c)

                gene_transcripts = gtf_obj.gene_trans_dt[gene_id]
                gene_coord_trans_dt[tuple(gene_c)].update(gene_transcripts)

            strand_B_gene_coordinates = sorted(strand_B_gene_coordinates)

            # Identify anti-sense transcripts in the 1st scaffold
            for mono_trans in strand_A_monoexonics:
                mono_c = gtf_obj.trans_exons_dt[mono_trans][0]
                for gene_c in strand_B_gene_coordinates:
                    # Important! This assumes gene coordinates are already sorted (gene_coord[0] < gene_coord[1])
                    if strict_nested:
                        if gene_c[0] <= mono_c[0] <= gene_c[1] and gene_c[0] <= mono_c[1] <= gene_c[1]:
                            monoexonic_antisense_set.add(mono_trans)

                            # Check the transcripts the model is antisense to
                            antisense_to = gene_coord_trans_dt[gene_c]
                            antisense_relation_dt[mono_trans].update(antisense_to)

                            # If the other model is also a monoexonic gene, then flag both as potential "antisense"
                            if len(antisense_to) == 1:
                                antisense_to = sorted(antisense_to)[0]
                                if len(gtf_obj.trans_exons_dt[antisense_to]) == 1:
                                    monoexonic_antisense_set.add(antisense_to)
                                    antisense_relation_dt[antisense_to].add(mono_trans)

                            if mono_c[0] > gene_c[0]:
                                # Using continue instead of break to capture all transcripts that are antisense with it
                                continue
                    else:
                        # Check for partial overlaps instead
                        if gene_c[0] <= mono_c[0] <= gene_c[1] or gene_c[0] <= mono_c[1] <= gene_c[1]:
                            monoexonic_antisense_set.add(mono_trans)

                            # Check the transcripts the model is antisense to
                            antisense_to = gene_coord_trans_dt[gene_c]
                            antisense_relation_dt[mono_trans].update(antisense_to)

                            # If the other model is also a monoexonic gene, then flag both as potential "antisense"
                            if len(antisense_to) == 1:
                                antisense_to = sorted(antisense_to)[0]
                                if len(gtf_obj.trans_exons_dt[antisense_to]) == 1:
                                    monoexonic_antisense_set.add(antisense_to)
                                    antisense_relation_dt[antisense_to].add(mono_trans)

                            if mono_c[0] > gene_c[0]:
                                continue

            # Identify anti-sense transcripts in the 2nd scaffold
            for mono_trans in strand_B_monoexonics:
                mono_c = gtf_obj.trans_exons_dt[mono_trans][0]
                for gene_c in strand_A_gene_coordinates:
                    # Important! This assumes gene coordinates are already sorted (gene_coord[0] < gene_coord[1])
                    if strict_nested:
                        if gene_c[0] <= mono_c[0] <= gene_c[1] and gene_c[0] <= mono_c[1] <= gene_c[1]:
                            monoexonic_antisense_set.add(mono_trans)

                            # Check the transcripts the model is antisense to
                            antisense_to = gene_coord_trans_dt[gene_c]
                            antisense_relation_dt[mono_trans].update(antisense_to)

                            # If the other model is also a monoexonic gene, then flag both as potential "antisense"
                            if len(antisense_to) == 1:
                                antisense_to = sorted(antisense_to)[0]
                                if len(gtf_obj.trans_exons_dt[antisense_to]) == 1:
                                    monoexonic_antisense_set.add(antisense_to)
                                    antisense_relation_dt[antisense_to].add(mono_trans)

                            if mono_c[0] > gene_c[0]:
                                continue
                                # break
                        else:
                            # Check for partial overlaps instead
                            if gene_c[0] <= mono_c[0] <= gene_c[1] or gene_c[0] <= mono_c[1] <= gene_c[1]:
                                monoexonic_antisense_set.add(mono_trans)

                                # Check the transcripts the model is antisense to
                                antisense_to = gene_coord_trans_dt[gene_c]
                                antisense_relation_dt[mono_trans].update(antisense_to)

                                # If the other model is also a monoexonic gene, then flag both as potential "antisense"
                                if len(antisense_to) == 1:
                                    antisense_to = sorted(antisense_to)[0]
                                    if len(gtf_obj.trans_exons_dt[antisense_to]) == 1:
                                        monoexonic_antisense_set.add(antisense_to)
                                        antisense_relation_dt[antisense_to].add(mono_trans)

                            if mono_c[0] > gene_c[0]:
                                continue
                                # break

    return monoexonic_antisense_set


def identify_antisense_fragments(gtf_obj, to_analyze=None, len_th=0.5):

    print(time.asctime(), f"Identifying antisense fragments")

    assert 0.0 <= len_th <= 1.0

    antisense_fragments, monoexon_fragments, intronic = [set() for _ in range(3)]

    if not to_analyze:
        to_analyze = set(gtf_obj.trans_exons_dt.keys())

    loci_group_dt = group_transcripts_across_strands(gtf_obj, to_analyze=to_analyze)

    # This dictionary is useful to write a table to explore the antisense fragments lengths
    len_ratio_dt = {}
    for gene_loci, overlap_group in loci_group_dt.items():
        longest_trans = get_longest_transcript(gtf_obj, overlap_group)
        longest_trans_strand = gtf_obj.trans_sense_dt[longest_trans]
        longest_len = get_transcript_length(gtf_obj, longest_trans)

        for t_id in sorted(overlap_group):
            # Analyze only potential antisense fragments in the overlap group, they are by definition monoexonic
            if len(gtf_obj.trans_exons_dt[t_id]) == 1:
                t_strand = gtf_obj.trans_sense_dt[t_id]

                t_len = get_transcript_length(gtf_obj, t_id)
                t_len_ratio = t_len / longest_len

                data = (t_len_ratio, t_len, longest_len)
                len_ratio_dt[t_id] = data

                # The comparison must be <= to remove those cases were the length is the same
                if t_len_ratio <= len_th and t_id != longest_trans:
                    if t_strand != longest_trans_strand:
                        antisense_fragments.add(t_id)
                    else:
                        monoexon_fragments.add(t_id)

        # Track intronic transcripts in case the user wants to keep them
        group_intronic = identify_intronic_transcripts(gtf_obj, overlap_group)
        intronic = intronic | group_intronic

    return antisense_fragments, monoexon_fragments, intronic, len_ratio_dt


def identify_unknown_chromosome_models(gtf_obj, genome_fa):

    print(time.asctime(), f'Identifying transcripts with unknown scaffold IDs', flush=True)

    print(time.asctime(), f'Extracting scaffold IDs from genome FASTA')
    genome_chrom = set()
    with open(genome_fa) as fh:
        for row in fh:
            if row.startswith(">"):
                # This allows to, most likely, get only the Chrom ID without any other information
                chrom_id = row.strip("\n").strip(">").replace("\t", "#").replace(" ", "#").split("#")[0]
                genome_chrom.add(chrom_id)

    unknown_chrom_transcripts = set()
    for chrom_id, chrom_trans in gtf_obj.chrom_trans_dt.items():
        # Important: GTF_OBJ Chrom IDs contain the strand information in the last position(Ex: chr01+)
        if chrom_id[:-1] not in genome_chrom:
            unknown_chrom_transcripts.update(chrom_trans)

    return unknown_chrom_transcripts


def remove_ambiguous_location_models(gtf_file, paths_dt, outname, to_add=None, to_keep=None, antisense_len=0, genome_fa=None):

    # print(time.asctime(), f'Analyzing transcripts models across strands for file: {gtf_file}', flush=True)

    if to_add is None:
        to_add = set()

    if to_keep is None:
        to_keep = set()

    gtf_obj = create_gtf_object(gtf_file)

    if genome_fa:
        unknown_chrom_transcripts = identify_unknown_chromosome_models(gtf_obj, genome_fa)
    else:
        unknown_chrom_transcripts = set()

    if 'unstranded' not in to_add:
        unstranded_transcripts = identify_unstranded_transcripts(gtf_obj)
    else:
        print(time.asctime(), 'Skipping identification of unstranded transcripts')
        unstranded_transcripts = set()

    monoexonic_antisense = identify_monoexonic_antisense(gtf_obj, to_ignore=unstranded_transcripts)
    antisense_fragments, monoexon_fragments, intronic, _ = identify_antisense_fragments(gtf_obj, monoexonic_antisense, antisense_len)

    if 'intronic' in to_add:
        monoexon_fragments = monoexon_fragments - intronic

    # Identify accepted transcripts
    tot_trans = set(gtf_obj.trans_exons_dt.keys())

    to_remove = unstranded_transcripts | antisense_fragments | monoexon_fragments | unknown_chrom_transcripts
    to_accept = tot_trans - to_remove

    outfile = write_gtf(gtf_obj, to_accept, paths_dt['inter'], f"{outname}_valid_loc.gtf")

    if 'additional' in to_keep:
        _ = write_gtf(gtf_obj, monoexonic_antisense, paths_dt['inter'], f"{outname}_antisense_monoexonic.gtf")

    if 'removed' in to_keep:
        if monoexonic_antisense:
            _ = write_gtf(gtf_obj, antisense_fragments, paths_dt['removed'], f"{outname}_fragments_antisense.gtf")
            _ = write_gtf(gtf_obj, monoexon_fragments, paths_dt['removed'], f"{outname}_fragments_monoexon.gtf")

        if unstranded_transcripts:
            _ = write_gtf(gtf_obj, unstranded_transcripts, paths_dt['removed'], f"{outname}_unstranded.gtf")

        if unknown_chrom_transcripts:
            _ = write_gtf(gtf_obj, unknown_chrom_transcripts, paths_dt['removed'], f"{outname}_unknown_chromosome.gtf")

    # Track numbers of accepted / removed
    global_track_dt[f"Removed#{outname}_fragments_antisense.gtf"] = len(antisense_fragments)
    global_track_dt[f"Removed#{outname}_fragments_monoexon.gtf"] = len(monoexon_fragments)
    global_track_dt[f"Removed#{outname}_unstranded.gtf"] = len(unstranded_transcripts)
    global_track_dt[f"Removed#{outname}_unknown_chromosome.gtf"] = len(unknown_chrom_transcripts)

    global_track_dt[f"Accepted#{outname}_valid_loc.gtf"] = len(to_accept)

    return outfile

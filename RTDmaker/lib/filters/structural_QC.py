import os
import time

from lib.tools.other_tools import *
from lib.tools.logger import logger
from lib.report.report_tools import global_track_dt
from lib.parsing.gtf_object_tools import create_gtf_object, write_gtf


def identify_groups_with_uniform_consensus(gtf_obj):

    # print(time.asctime(), "Identifying genomic regions with a single-consensus group")

    global_single_consensus_group = set()
    for chrom, strand_transcripts in gtf_obj.chrom_trans_dt.items():

        # Group transcripts by their overlap
        overlapping_transcripts = group_transcripts_by_overlap(gtf_obj, strand_transcripts)

        # Group transcripts by their boundary consensus (if their first/last exon overlap)
        consensus_groups_list, single_transcripts = [], set()
        for overlap_group in overlapping_transcripts:
            if len(overlap_group) < 2:
                single_transcripts.update(overlap_group)
            else:
                consensus_groups = group_transcripts_by_consensus(gtf_obj, overlap_group)
                consensus_groups_list.append(consensus_groups)

        # Genes containing a single transcript
        global_single_consensus_group.update(single_transcripts)

        # Identify if the overlapping region have multiple consensus groups
        multiple_consensus_groups, one_consensus_group = [], set()
        for consensus_groups in consensus_groups_list:
            if len(consensus_groups) == 1:
                one_consensus_group.update(flat(consensus_groups))
            else:
                multiple_consensus_groups.append(consensus_groups)

        # Overlapping groups containing a single boundary consensus group
        global_single_consensus_group.update(one_consensus_group)

    return global_single_consensus_group


def identify_non_overlapping_subgroups(gtf_obj, to_analyze):

    # print(time.asctime(), "Identifying genomic regions with multiple boundary-consensus subgroups")

    group_with_segments, group_without_segments = [set() for _ in range(2)]
    for chrom, strand_transcripts in gtf_obj.chrom_trans_dt.items():

        strand_transcripts = [t_id for t_id in strand_transcripts if t_id in to_analyze]

        # Group transcripts by their overlap
        overlapping_transcripts = group_transcripts_by_overlap(gtf_obj, strand_transcripts)

        # Group transcripts by their boundary consensus (if their first/last exon overlap)
        for overlap_group in overlapping_transcripts:
            consensus_groups = group_transcripts_by_consensus(gtf_obj, overlap_group)

            all_subgroups_overlap = True
            # group_with_segments, group_without_segments = (set() for _ in range(2))
            for (group_A, group_B) in combinations(consensus_groups, 2):

                trans_A = get_longest_transcript(gtf_obj, group_A)
                trans_B = get_longest_transcript(gtf_obj, group_B)

                trans_A_boundaries = get_start_end(get_exon_model(gtf_obj, trans_A))
                trans_B_boundaries = get_start_end(get_exon_model(gtf_obj, trans_B))

                # Check if the overlap group contain two models that don't overlap between them ("segmented groups")
                if not get_overlap_percentage(trans_A_boundaries, trans_B_boundaries):
                    all_subgroups_overlap = False

            if all_subgroups_overlap:
                for group in consensus_groups:
                    group_without_segments.update(group)
            else:
                for group in consensus_groups:
                    group_with_segments.update(group)

    return group_with_segments, group_without_segments


def identify_intronic_transcripts(gtf_obj, t_group):

    trans_exons_dt = gtf_obj.trans_exons_dt

    intronic_trans = set()

    if len(t_group) < 2:
        return intronic_trans

    # Get all the exon-ranges in the group and identify which transcripts are mono-exonic
    monoexons, group_exons = [], set()
    for t_id in t_group:
        trans_exons = trans_exons_dt[t_id]

        if len(trans_exons) < 2:
            monoexons.append(t_id)
        else:
            group_exons.update(trans_exons)

    group_exons = sorted(group_exons, key=lambda exon: exon[0])

    for mono_trans in monoexons:
        overlap_exon = False
        mono_trans_model = get_exon_model(gtf_obj, mono_trans)

        for exon in group_exons:
            if get_overlap_percentage(mono_trans_model[0], exon):
                overlap_exon = True
                break

        if not overlap_exon:
            intronic_trans.add(mono_trans)

    return intronic_trans


def identify_overlapping_transcripts(gtf_obj, overlap_group):

    accepted_overlapping, rejected_overlapping, rejected_segmented = [set() for _ in range(3)]

    # Overlapping identification consist of the following checks:
    # 1) Identify pair of "consensus groups" that do NOT overlap
    # 2) Identify transcripts that overlap the pair of "non-overlapping consensus groups"
    # 3) Check if the number of overlapping transcripts is more than the transcripts that they overlap
    # 4) If so, accept the overlapping transcripts (and ID the rest as "segments"), otherwise remove the overlap trans

    # Identify consensus groups
    consensus_groups = group_transcripts_by_consensus(gtf_obj, overlap_group)

    # Ignore "mono-exonic isoforms" when performing this comparison as they are more likely to be assembler artifacts
    valid_consensus_groups = []
    for group in consensus_groups:
        if len(group) == 1:
            trans_id = group[0]
            trans_exons = gtf_obj.trans_exons_dt[trans_id]
            if len(trans_exons) == 1:
                # Ignore these "single mono-exonic isoforms"
                continue
        valid_consensus_groups.append(group)

    # Identify pair of non-overlapping consensus groups
    non_overlapping_group_pairs = []
    for (group_A, group_B) in combinations(valid_consensus_groups, 2):
        # Get intron and exon boundaries
        group_A_exon_boundaries, group_A_intron_boundaries = get_group_boundaries(gtf_obj, group_A)
        group_B_exon_boundaries, group_B_intron_boundaries = get_group_boundaries(gtf_obj, group_B)

        # Check if the two groups DON'T overlap
        if not get_overlap_percentage(group_A_exon_boundaries, group_B_intron_boundaries):
            if not get_overlap_percentage(group_A_intron_boundaries, group_B_exon_boundaries):
                group_pair = (group_A, group_B)
                non_overlapping_group_pairs.append(group_pair)

    overlapping_groups_dt = defaultdict(set)
    potential_overlapping = set()
    for (group_A, group_B) in non_overlapping_group_pairs:
        group_A_exon_boundaries, group_A_intron_boundaries = get_group_boundaries(gtf_obj, group_A)
        group_B_exon_boundaries, group_B_intron_boundaries = get_group_boundaries(gtf_obj, group_B)

        # Check if there is a 3rd group that overlap both not-overlapping groups
        for test_group in valid_consensus_groups:
            if test_group != group_A and test_group != group_B:
                test_exon_boundaries, test_intron_boundaries = get_group_boundaries(gtf_obj, test_group)

                if get_overlap_percentage(test_exon_boundaries, group_A_exon_boundaries):
                    if get_overlap_percentage(test_exon_boundaries, group_B_exon_boundaries):
                        potential_overlapping.add(tuple(test_group))

                        # # Keep track of the group that overlap with it
                        # overlapping_groups_dt[tuple(test_group)].add(tuple(group_A))
                        # overlapping_groups_dt[tuple(test_group)].add(tuple(group_B))

                        # Keep track of the group that overlap with it
                        overlap_pair = (tuple(group_A), tuple(group_B))
                        overlapping_groups_dt[tuple(test_group)].add(overlap_pair)

    # Check if the number of overlapping transcripts in a group are more than the number of transcripts they overlap
    for over_group in sorted(potential_overlapping):

        # Transcripts in the overlapping group
        overlapping_trans = [t_id for t_id in over_group]

        # Transcripts of the other groups that overlap with the group being analyze
        # over_group is already a tuple, no need to convert it to use as key
        overlap_subgroups = overlapping_groups_dt[over_group]

        for (group_A, group_B) in overlap_subgroups:

            subgroups_trans = []
            for subgroup in (group_A, group_B):
                for t_id in subgroup:
                    if t_id not in subgroups_trans:
                        subgroups_trans.append(t_id)

            # Adding +1 to the subgroups favor the overlapping transcripts over the segments when the numbers are equal
            if len(overlapping_trans) + 1 >= len(subgroups_trans):
                accepted_overlapping.update(overlapping_trans)
                rejected_segmented.update(subgroups_trans)
            else:
                # If there is a pair of non-overlapping subgroups more numerous, reject the current overlapping models
                rejected_overlapping.update(overlapping_trans)

                # Re-assign the previously classified models
                for t_id in overlapping_trans:
                    try:
                        accepted_overlapping.remove(t_id)
                    except KeyError:
                        continue

                for t_id in subgroups_trans:
                    try:
                        rejected_overlapping.remove(t_id)
                    except KeyError:
                        continue
                # If the overlapping transcripts fails once, there is no need to test other cases
                break

    # Ensure that the classifications are unique
    accepted_overlapping = accepted_overlapping - rejected_overlapping

    # Identify and recover intronic transcripts (We don't want to count them as chimeric/fragments)
    intronic = identify_intronic_transcripts(gtf_obj, overlap_group)

    rejected_segmented = rejected_segmented - intronic

    return accepted_overlapping, rejected_overlapping, rejected_segmented


def identify_fragmentary_transcripts(gtf_obj, overlap_group, len_th):

    fragments, overextended = [set() for _ in range(2)]

    # Identify and recover intronic transcripts to avoid classifying them as fragments
    intronic = identify_intronic_transcripts(gtf_obj, overlap_group)

    # Get longest transcrpipt
    longest_len = get_longest_length(gtf_obj, overlap_group)

    potential_fragments, potential_accepted = [set() for _ in range(2)]
    for trans_id in overlap_group:
        trans_len = get_transcript_length(gtf_obj, trans_id)
        # Calculate length ratio
        if (trans_len / longest_len) >= len_th:
            potential_accepted.add(trans_id)
        else:
            potential_fragments.add(trans_id)

    # Check for overextended transcripts. Overextend transcrips is defined by:
    # a) the accepted transcript is only 1,
    # b) and there is at least 1 consensus subgroup with multiple transcripts
    fragment_consensus_groups = group_transcripts_by_consensus(gtf_obj, potential_fragments)
    if len(potential_accepted) == 1:
        for fragment_subgroup in fragment_consensus_groups:
            if len(fragment_subgroup) >= len(potential_accepted) * 2:
                overextended.update(potential_accepted)

    if not overextended:
        fragments.update(potential_fragments)
        fragments = fragments - intronic
        return fragments, overextended

    # If there is an overextended transcript, perform a new length check without the overextended transcript
    subgroup = [t_id for t_id in overlap_group if t_id not in overextended]

    # Re-initialize the fragment group
    fragments = set()

    # Keep the longest overlapping group, flag other shorter non-overlapping subgroups as fragments
    # Beware! The sorted(reverse=True) is giving incorrect results, do NOT use it!
    # Why is this happening? Because ¯\_(ツ)_/¯
    sorted_subgroups = sorted(group_transcripts_by_overlap(gtf_obj, subgroup),
                              key=lambda group: get_longest_length(gtf_obj, group))

    selected_subgroup = sorted_subgroups[0]
    fragments.update(flat(sorted_subgroups[:-1]))

    ref_len = get_longest_length(gtf_obj, selected_subgroup)
    for trans_id in selected_subgroup:
        trans_len = get_transcript_length(gtf_obj, trans_id)

        # We don't want to re-introduce fragments that are overly short; perform final, more tolerant, length check
        soft_len_th = 0.33
        if (trans_len / longest_len) >= soft_len_th:
            if (trans_len / ref_len) < len_th:
                fragments.add(trans_id)
        else:
            fragments.add(trans_id)

    fragments = fragments - intronic

    return fragments, overextended


def overlap_analysis(gtf_obj, to_analyze=None, write_table=False, verb=False):

    if verb:
        print(time.asctime(), "Identifying overlapping transcripts", flush=True)

    if write_table:
        fh = open("chimeric_table.csv", "w+")
        fh.write("Potential_chimeric,Potential_fragments\n")

    global_accepted_overlapping, global_rejected_overlapping, global_rejected_segments = [set() for _ in range(3)]
    for chrom, strand_transcripts in gtf_obj.chrom_trans_dt.items():

        # Filter out transcripts from groups that don't contain overlapping transcripts
        if to_analyze:
            strand_transcripts = [t_id for t_id in strand_transcripts if t_id in to_analyze]

        # Group transcripts by their overlap
        overlapping_transcripts = group_transcripts_by_overlap(gtf_obj, strand_transcripts)

        # Group transcripts by their boundary consensus (if their first/last exon overlap)
        for overlap_group in overlapping_transcripts:
            accepted_over, rejected_over, rejected_seg = identify_overlapping_transcripts(gtf_obj, overlap_group)

            if write_table and rejected_over:
                    pot_chimeric = ";".join(sorted(rejected_over))
                    pot_fragment = ";".join(sorted(accepted_over))
                    fh.write(f"{pot_chimeric},{pot_fragment}\n")

            # Update clasifications
            global_accepted_overlapping = global_accepted_overlapping | accepted_over
            global_rejected_overlapping = global_rejected_overlapping | rejected_over
            global_rejected_segments = global_rejected_segments | rejected_seg

    if write_table:
        fh.close()

    return global_accepted_overlapping, global_rejected_overlapping, global_rejected_segments


def fragmentary_analysis(gtf_obj, to_analyze, len_th, verb=False):

    if verb:
        print(time.asctime(), "Identifying fragmentary transcripts", flush=True)

    global_fragments, global_overextended = [set() for _ in range(2)]
    for chrom, strand_transcripts in gtf_obj.chrom_trans_dt.items():

        # Filter out transcripts to ignore when checking for fragments, Example:
        # Transcripts from groups without fragments, and rejected overlapping transcripts models
        strand_transcripts = [t_id for t_id in strand_transcripts if t_id in to_analyze]

        # Group transcripts by their overlap
        overlapping_transcripts = group_transcripts_by_overlap(gtf_obj, strand_transcripts)

        # Group transcripts by their boundary consensus (if their first/last exon overlap)
        for overlap_group in overlapping_transcripts:

            fragments, overextended = identify_fragmentary_transcripts(gtf_obj, overlap_group, len_th)

            global_fragments = global_fragments | fragments
            global_overextended = global_overextended | overextended

    return global_fragments, global_overextended


def remove_chimeric_and_fragments(gtf_file, len_th, paths_dt, outname, logfile, to_add=None, to_keep=None):

    logger(logfile)

    print(time.asctime(), f"Removing chimeric and fragmentary transcripts", flush=True)

    if to_add is None:
        to_add = set()

    if to_keep is None:
        to_keep = set()

    gtf_obj = create_gtf_object(gtf_file)

    # All transcripts in the annotation
    all_set = set(gtf_obj.trans_exons_dt.keys())
    fragments = set()

    print(time.asctime(), f"Classifying transcripts into groups according to their exonic boundaries", flush=True)
    # Identifying regions with a single consensus groups
    uniform_consensus = identify_groups_with_uniform_consensus(gtf_obj)

    # Identify groups where there are multiple subgroups of non-overlapping boundary-consensus groups
    contain_subgroups = all_set - uniform_consensus
    contain_overlaps, no_overlaps = identify_non_overlapping_subgroups(gtf_obj, contain_subgroups)

    print(time.asctime(), f"Identifying overlapping transcripts", flush=True)
    # Identify overlapping transcripts to allow for an accurate identification of fragmentary transcripts
    potential_accepted, potential_overlap, potential_segment = overlap_analysis(gtf_obj, contain_overlaps)

    print(time.asctime(), f"Removing fragmentary transcripts", flush=True)

    # Groups to ignore while checking for fragments,
    contain_fragments = all_set - potential_overlap

    temp_fragments, overextended = fragmentary_analysis(gtf_obj, contain_fragments, len_th)
    fragments = fragments | temp_fragments

    # It is possible to perform a 2nd analysis iteration
    refine_models = True
    if refine_models:
        # If it contains subgroups, but no overlaps, check only for fragments
        temp_fragments, _ = fragmentary_analysis(gtf_obj, no_overlaps, len_th)
        fragments = fragments | temp_fragments

        # Remove fragments from potentially accepted overlap groups
        temp_fragments, _ = fragmentary_analysis(gtf_obj, potential_accepted, len_th)
        fragments = fragments | temp_fragments

        # Remove fragments from the potentially chimeric models, this make their visualization on IGB/IGV more clean
        temp_fragments, _ = fragmentary_analysis(gtf_obj, potential_overlap, len_th)
        fragments = fragments | temp_fragments

        temp_fragments, _ = fragmentary_analysis(gtf_obj, potential_segment, len_th)
        fragments = fragments | temp_fragments

        # Finally, perform the overlap analysis to decide whether to keep or remove the overlapping transcripts
        no_fragments = contain_overlaps - fragments
        selected_overlap, rejected_overlap, rejected_segment = overlap_analysis(gtf_obj, no_fragments)

        # Group the fragmentary transcripts together to annotate them into the same output file
        to_check = all_set - (uniform_consensus | rejected_overlap | rejected_segment | fragments)

        temp_fragments, _ = fragmentary_analysis(gtf_obj, to_check, len_th)
        fragments = fragments | temp_fragments

    else:
        selected_overlap, rejected_overlap, rejected_segment = [set() for _ in range(3)]
        rejected_overlap = rejected_overlap | potential_overlap
        rejected_segment = rejected_segment | potential_segment

    # Refine the transcript classifications to make output easier to explore with IGB/IGV
    refined_overlaps = potential_accepted - (rejected_overlap | rejected_segment | fragments)
    accepted_overlap = (selected_overlap | refined_overlaps)

    # Identify the subset of accepted overlaps that were recovered
    selected_overlap = selected_overlap - (potential_overlap | fragments)

    # Remove potential fragments from the other removed groups to avoid information redundancy
    rejected_overlap = rejected_overlap - (rejected_segment | fragments)
    rejected_segment = rejected_segment - (rejected_overlap | fragments)

    # Allow the user to remove intronic transcripts if they wish
    intronic = set()
    if 'intronic' not in to_add:
        print(time.asctime(), 'Removing intronic transcripts')
        for chrom, strand_transcripts in gtf_obj.chrom_trans_dt.items():
            overlapping_transcripts = group_transcripts_by_overlap(gtf_obj, strand_transcripts)
            for overlap_group in overlapping_transcripts:
                intronic_group = identify_intronic_transcripts(gtf_obj, overlap_group)
                intronic = intronic | intronic_group

    fragments = fragments | intronic

    # Remove the rejected transcript models
    to_remove = fragments | rejected_overlap | rejected_segment
    accepted = set(gtf_obj.trans_gene_dt.keys()) - to_remove

    outfile = write_gtf(gtf_obj, accepted, paths_dt["inter"], f"{outname}_to_reannotate.gtf")

    if 'additional' in to_keep:
        _ = write_gtf(gtf_obj, overextended, paths_dt["additional"], f"{outname}_overextended.gtf")
        _ = write_gtf(gtf_obj, accepted_overlap, paths_dt["additional"], f"{outname}_accepted_overlap.gtf")
        _ = write_gtf(gtf_obj, rejected_overlap, paths_dt["additional"], f"{outname}_rejected_overlap.gtf")

        if selected_overlap:
            _ = write_gtf(gtf_obj, selected_overlap, paths_dt["additional"], f"{outname}_refined_overlap.gtf")

    if 'removed' in to_keep:
        _ = write_gtf(gtf_obj, rejected_overlap, paths_dt["removed"], f"{outname}_rejected_overlap.gtf")
        _ = write_gtf(gtf_obj, rejected_segment, paths_dt["removed"], f"{outname}_rejected_segment.gtf")
        _ = write_gtf(gtf_obj, fragments, paths_dt["removed"], f"{outname}_fragments.gtf", w_mode="a+")

    # Track numbers of accepted / removed
    global_track_dt[f"Removed#{outname}_rejected_overlap.gtf"] = len(rejected_overlap)
    global_track_dt[f"Removed#{outname}_rejected_segment.gtf"] = len(rejected_segment)
    global_track_dt[f"Removed#{outname}_fragments.gtf"] = len(fragments)
    global_track_dt[f"Accepted#{outname}_to_reannotate.gtf"] = len(accepted)

    return outfile

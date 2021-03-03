import os
import time
import warnings
from collections import defaultdict

from lib.filters.structural_QC import overlap_analysis
from lib.filters.structural_QC import identify_groups_with_uniform_consensus, identify_non_overlapping_subgroups

from lib.filters.GenLoc_QC import identify_monoexonic_antisense, identify_antisense_fragments, get_longest_transcript
from lib.filters.structural_QC import group_transcripts_by_overlap, identify_intronic_transcripts, group_transcripts_by_consensus

from lib.tools.logger import logger
from lib.parsing.gtf_object_tools import create_gtf_object

from lib.report.report_tools import percentiles_from_counts, simple_write_table, get_per, get_distribution_count

warnings.filterwarnings("ignore")


def generate_annotation_report(gtf_file, outpath, outname, logfile):

    # TODO handle error if the input file is not a GTF file

    logger(logfile)

    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    print(time.asctime(), f"Generating report of {gtf_file}", flush=True)

    gtf_obj = create_gtf_object(gtf_file)

    # Analyze transcriptome annotation (Number of genes/transcripts, etc)
    report_folder = gtf_exploration(gtf_obj, outpath, outname)

    print(time.asctime(), f"Report tables location: {report_folder}", flush=True)


def report_chimeric(gtf_obj, outpath, outname):

    print(time.asctime(), f"Identifying chimeric transcripts", flush=True)
    # All transcripts in the annotation
    all_set = set(gtf_obj.trans_exons_dt.keys())

    # print(time.asctime(), f"Classifying transcripts into groups according to their exonic boundaries", flush=True)
    # Identifying regions with a single consensus groups
    uniform_consensus = identify_groups_with_uniform_consensus(gtf_obj)

    # Identify groups where there are multiple subgroups of non-overlapping boundary-consensus groups
    contain_subgroups = all_set - uniform_consensus
    contain_overlaps, no_overlaps = identify_non_overlapping_subgroups(gtf_obj, contain_subgroups)

    # print(time.asctime(), f"Identifying overlapping transcripts", flush=True)
    # Identify overlapping transcripts to allow for an accurate identification of fragmentary transcripts
    potential_accepted, potential_overlap, potential_segment = overlap_analysis(gtf_obj, contain_overlaps)

    # Generate rows for table
    # First relate chimeric trans with their genes
    potential_chimeric = potential_accepted | potential_overlap

    # Chimeric genes/transcripts
    chimeric_genes_dt = defaultdict(set)
    for t_id in potential_chimeric:
        g_id = gtf_obj.trans_gene_dt[t_id]
        chimeric_genes_dt[g_id].add(t_id)

    chimeric_table_rows = ["Gene ID\tChimeric transcripts\n"]
    for g_id, g_trans in sorted(chimeric_genes_dt.items(), key=lambda kv: len(kv[1]), reverse=True):
        trans_str = ";".join(sorted(g_trans))
        row = f"{g_id}\t{trans_str}\n"
        chimeric_table_rows.append(row)

    # Split genes/transcripts
    split_genes_dt = defaultdict(set)
    for t_id in potential_segment:
        g_id = gtf_obj.trans_gene_dt[t_id]
        split_genes_dt[g_id].add(t_id)

    split_table_rows = ["Gene ID\tSplit Transcripts\n"]
    for g_id, g_trans in sorted(split_genes_dt.items(), key=lambda kv: len(kv[1]), reverse=True):
        trans_str = ";".join(sorted(g_trans))
        row = f"{g_id}\t{trans_str}\n"
        split_table_rows.append(row)

    # Write tables
    table_data = [("chimeric", chimeric_table_rows), ("split", split_table_rows)]
    for cat_tag, cat_rows in table_data:
        outfile = os.path.join(outpath, f"{outname}_potentially_{cat_tag}_genes.tsv")
        simple_write_table(cat_rows, outfile)


def report_not_uniform_groups(gtf_obj, outpath, outname, to_ignore=None, ratio_th=0.7):

    not_uniform_transcripts = set()

    print(time.asctime(), f"Identifying genes potentially containing fragments", flush=True)
    # All transcripts in the annotation
    all_set = set(gtf_obj.trans_exons_dt.keys())

    if not to_ignore:
        to_ignore = set()

    all_set = all_set - to_ignore

    # print(time.asctime(), f"Classifying transcripts into groups according to their exonic boundaries", flush=True)
    # Identifying regions with a single consensus groups
    uniform_consensus = identify_groups_with_uniform_consensus(gtf_obj)

    # Identify groups where there are multiple subgroups
    contain_subgroups = all_set - uniform_consensus

    # contain_overlaps, no_overlaps = identify_non_overlapping_subgroups
    flat = lambda l: [e for sub in l for e in sub]

    for chrom, strand_transcripts in gtf_obj.chrom_trans_dt.items():

        # Transcripts to analyze
        strand_transcripts = [t_id for t_id in strand_transcripts if t_id in contain_subgroups]

        # Remove mono-exonic transcripts from the analysis to avoid noise (due to intronic, redundant monoexons, etc)
        strand_transcripts = [t_id for t_id in strand_transcripts if len(gtf_obj.trans_exons_dt[t_id]) > 1]

        # Group transcripts by their overlap
        overlapping_transcripts = group_transcripts_by_overlap(gtf_obj, strand_transcripts)

        # Group transcripts by their boundary consensus (if their first/last exon overlap)
        for overlap_group in overlapping_transcripts:
            consensus_groups = group_transcripts_by_consensus(gtf_obj, overlap_group)

            n_groups = len(consensus_groups)

            group_trans = flat(consensus_groups)
            n_trans = len(group_trans)

            group_to_trans_ratio = n_groups / n_trans

            # This ratio is high if most transcripts in the group are in their own boundary-consensus subgroup
            if group_to_trans_ratio >= ratio_th:
                # Ignore groups with less than 3 isoforms (no sense reporting genes with only 2 isoforms with diff ends)
                if n_trans > 2:
                    not_uniform_transcripts.update(group_trans)

    # Chimeric genes/transcripts
    genes_trans_dt = defaultdict(set)
    for t_id in not_uniform_transcripts:
        g_id = gtf_obj.trans_gene_dt[t_id]
        genes_trans_dt[g_id].add(t_id)

    not_uniform_table_rows = ["Gene ID\tTranscripts IDs\n"]
    for g_id, g_trans in sorted(genes_trans_dt.items(), key=lambda kv: len(kv[1]), reverse=True):
        trans_str = ";".join(sorted(g_trans))
        row = f"{g_id}\t{trans_str}\n"
        not_uniform_table_rows.append(row)

    not_uniform_table = os.path.join(outpath, f"{outname}_potentially_fragmented_genes.tsv")
    simple_write_table(not_uniform_table_rows, not_uniform_table)


def gtf_exploration(gtf_obj, outpath, outname):

    # Calculating Nº of Scaffolds/Genes/Transcripts
    outfile = os.path.join(outpath, f"{outname}")

    # 0) Get number of genes and transcripts to calculate percentages
    n_chrom = len(gtf_obj.chrom_gene_dt)
    n_genes = len(gtf_obj.gene_coords_dt)
    n_trans = len(gtf_obj.trans_exons_dt)
    p_genes = get_per(n_genes, n_genes)
    p_trans = get_per(n_trans, n_trans)

    # 1) Scaffold table
    header = "Scaffold\tNº of Genes\tGenes %\tNº of Transcripts\tTranscripts %\n"
    scaffold_rows = [header]

    strand_genes_dt, strand_trans_dt = [defaultdict(set) for _ in range(2)]
    chrom_genes_dt, chrom_trans_dt = [{} for _ in range(2)]
    for chrom in gtf_obj.chrom_gene_dt.keys():
        n_chrom_genes = len(gtf_obj.chrom_gene_dt[chrom])
        chrom_genes_dt[chrom] = n_chrom_genes

        n_chrom_trans = len(gtf_obj.chrom_trans_dt[chrom])
        chrom_trans_dt[chrom] = n_chrom_trans

        chrom_strand = chrom[-1]
        strand_genes_dt[chrom_strand].update(gtf_obj.chrom_gene_dt[chrom])
        strand_trans_dt[chrom_strand].update(gtf_obj.chrom_trans_dt[chrom])

    for chrom in sorted(chrom_genes_dt.keys()):
        n_chrgenes = chrom_genes_dt[chrom]
        p_chrgenes = get_per(n_chrgenes, n_genes)

        n_chrtrans = chrom_trans_dt[chrom]
        p_chrtrans = get_per(n_chrtrans, n_trans)

        row = f"{chrom}\t{n_chrgenes:>7}\t{p_chrgenes}\t{n_chrtrans:>7}\t{p_chrtrans}\n"
        scaffold_rows.append(row)

    unstranded_genes = set()
    for strand, genes_set in strand_genes_dt.items():
        n_strand_genes = len(genes_set)
        p_strand_genes = get_per(n_strand_genes, n_genes)

        n_strand_trans = len(strand_trans_dt[strand])
        p_strand_trans = get_per(n_strand_trans, n_trans)

        # Ambiguous strand
        if strand not in {"+", "-"}:
            unstranded_genes.update(genes_set)

        row = f"Strand {strand}\t{n_strand_genes:>7}\t{p_strand_genes}\t{n_strand_trans:>7}\t{p_strand_trans}\n"
        scaffold_rows.append(row)

    # Write table
    simple_write_table(scaffold_rows, f"{outfile}_numbers_scaffolds.tsv")

    # 2) Gene category table
    gene_categories_rows = []

    header = "Category\tNº of Genes\tGenes %\tNº of Transcripts\tTranscripts %\n"
    totals = f"Total\t{n_genes}\t{p_genes}\t{n_trans}\t{p_trans}\n"
    gene_categories_rows.extend([header, totals])

    monoexon_genes, single_multiexon_genes, multi_iso_genes = [set() for _ in range(3)]
    monoexon_trans, single_multiexon_trans, multi_iso_trans = [set() for _ in range(3)]

    n_isoform_dt = defaultdict(set)
    categories_data = []
    for gene_id, gene_transcripts in gtf_obj.gene_trans_dt.items():
        if len(gene_transcripts) < 2:
            t_id = sorted(gene_transcripts)[0]
            trans_exons = gtf_obj.trans_exons_dt[t_id]
            if len(trans_exons) < 2:
                monoexon_genes.add(gene_id)
                monoexon_trans.update(gene_transcripts)

            else:
                single_multiexon_genes.add(gene_id)
                single_multiexon_trans.update(gene_transcripts)

        else:
            multi_iso_genes.add(gene_id)
            multi_iso_trans.update(gene_transcripts)

        len_key = len(gene_transcripts)
        n_isoform_dt[len_key].add(gene_id)

    categories_data.append(("Multi-Isoforms", multi_iso_genes, multi_iso_trans))
    categories_data.append(("Single-Isoform Intron-Containing", single_multiexon_genes, single_multiexon_trans))
    categories_data.append(("Monoexonic", monoexon_genes, monoexon_trans))

    for cat_name, cat_genes, cat_trans in categories_data:

        n_cat_genes = len(cat_genes)
        p_cat_genes = get_per(n_cat_genes, n_genes)

        n_cat_trans = len(cat_trans)
        p_cat_trans = get_per(n_cat_trans, n_trans)

        gene_categories_rows.append(f"{cat_name}\t{n_cat_genes}\t{p_cat_genes}\t{n_cat_trans}\t{p_cat_trans}\n")

    # Write table
    simple_write_table(gene_categories_rows, f"{outfile}_numbers_gene_categories.tsv")

    # 2.5) Find additional gene categories of interest

    other_categories_data = []

    # Unstranded genes
    unstranded_trans = set()
    for g_id in unstranded_genes:
        unstranded_trans.update(gtf_obj.gene_trans_dt[g_id])

    other_categories_data.append(("Unstranded", unstranded_genes, unstranded_trans))

    # Identify Intronic and Redundant Mono-exonic transcripts
    intronic_genes, intronic_trans, red_mono_genes, red_mono_trans = [set() for _ in range(4)]
    for chrom, strand_transcripts in gtf_obj.chrom_trans_dt.items():
        overlapping_transcripts = group_transcripts_by_overlap(gtf_obj, strand_transcripts)
        for overlap_group in overlapping_transcripts:
            # In the case of references, this group contains both intronic and redundant monoexonic transcripts
            temp_group = identify_intronic_transcripts(gtf_obj, overlap_group)
            t_group = set(temp_group)

            for t_id in sorted(t_group):
                g_id = gtf_obj.trans_gene_dt[t_id]

                # If all the transcripts in the gene are redundant, then tag it as a redundant monoexonic gene
                g_trans = gtf_obj.gene_trans_dt[g_id]
                if set(g_trans).difference(t_group):
                    intronic_genes.add(g_id)
                    intronic_trans.add(t_id)
                else:
                    t_longest = get_longest_transcript(gtf_obj, t_group)
                    if t_id != t_longest:
                        red_mono_trans.add(t_id)
                        red_mono_genes.add(g_id)

    other_categories_data.append(("Intronic", intronic_genes, intronic_trans))
    other_categories_data.append(("Monoexonic Redundant", red_mono_genes, red_mono_trans))

    # Identify Monoexonic Antisense Genes
    # Set length value to 100% so it identifies all antisense transcripts
    antisense = identify_monoexonic_antisense(gtf_obj)
    antisense_trans, _, _, _ = identify_antisense_fragments(gtf_obj, to_analyze=antisense, len_th=1.0)

    antisense_genes = set()
    for t_id in antisense_trans:
        g_id = gtf_obj.trans_gene_dt[t_id]
        antisense_genes.add(g_id)

    other_categories_data.append(("Monoexonic Antisense", antisense_genes, antisense_trans))

    header = "Category\tNº of Genes\tGenes %\tNº of Transcripts\tTranscripts %\n"
    other_categories_rows = [header]
    for cat_name, cat_genes, cat_trans in other_categories_data:

        n_cat_genes = len(cat_genes)
        p_cat_genes = get_per(n_cat_genes, n_genes)

        n_cat_trans = len(cat_trans)
        p_cat_trans = get_per(n_cat_trans, n_trans)

        other_categories_rows.append(f"{cat_name}\t{n_cat_genes}\t{p_cat_genes}\t{n_cat_trans}\t{p_cat_trans}\n")

    # Write table
    simple_write_table(other_categories_rows, f"{outfile}_numbers_other_categories.tsv")

    # Write Transcript IDs of other categories into a table
    for cat_name, cat_genes, cat_trans in other_categories_data:
        table_name = cat_name.lower().replace(" ", "_").replace("-", "")
        temp_rows = ["Transcript_ID\n"]
        for t_id in sorted(cat_trans):
            temp_rows.append(f"{t_id}\n")

        if temp_rows:
            simple_write_table(temp_rows, f"{outfile}_IDs_other_categories_{table_name}.tsv")

    # 3) Number of Isoforms table
    iso_table_rows = []

    header = "Nº of Isoforms\tNº of Genes\tGenes %\tNº of Transcripts\tTranscripts %\n"
    iso_table_rows.append(header)

    iso_id_rows = ["Gene ID\tNº of Isoforms\n"]

    flat = lambda l: [e for sub in l for e in sub]

    # TODO add a +20/+n category?
    for n_iso, gene_set in sorted(n_isoform_dt.items()):

        trans_set = set(flat([gtf_obj.gene_trans_dt[g_id] for g_id in gene_set]))

        n_iso_genes = len(gene_set)
        p_iso_genes = get_per(n_iso_genes, n_genes)

        n_iso_trans = len(trans_set)
        p_iso_trans = get_per(n_iso_trans, n_trans)

        row = f"{n_iso}\t{n_iso_genes}\t{p_iso_genes}\t{n_iso_trans}\t{p_iso_trans}\n"
        iso_table_rows.append(row)

        for g_id in sorted(gene_set):
            id_row = f"{g_id}\t{n_iso}\n"
            iso_id_rows.append(id_row)

    # Write tables
    simple_write_table(iso_table_rows, f"{outfile}_numbers_isoforms.tsv")
    simple_write_table(iso_id_rows, f"{outfile}_IDs_isoforms.tsv")

    # 3.5) Report potential chimeric present in the annotation
    report_chimeric(gtf_obj, outpath, outname)

    # 3.75) Report genes and transcripts with large number of BC subgroups (may indicate fragmentation)
    report_not_uniform_groups(gtf_obj, outpath, outname, to_ignore=intronic_trans | red_mono_trans)

    # 4) Lengths tables
    # Dict for the lengths counts to use a more memory efficient approach to generate percentiles
    trans_len_count, exons_len_count, introns_len_count, cds_len_count, monoexon_len_count = [{} for _ in range(5)]

    # Dict to store the relationship between Transcript ID and exon/introns to generate table downstream
    trans_lenght_dt, trans_exons_dt, trans_introns_dt = [defaultdict(set) for _ in range(3)]

    # Temporary variables to store the lengths distribution to generate dict counts for these categories
    trans_lengths, cds_lengths, monoexon_lengths = [[] for _ in range(3)]

    for trans_id, trans_exons in gtf_obj.trans_exons_dt.items():
        if len(trans_exons) > 1:
            trans_introns = gtf_obj.trans_introns_dt[trans_id]
        else:
            trans_introns = []

        # Transcripts lengths
        trans_len = (trans_exons[-1][-1] - trans_exons[0][0]) + 1
        trans_lengths.append(trans_len)

        # Mono-exon genes lengths
        if trans_id in monoexon_trans:
            monoexon_lengths.append(trans_len)

        trans_lenght_dt[trans_id] = trans_len

        # CDS lengths
        trans_cds = gtf_obj.trans_cds_dt[trans_id]
        if trans_cds:
            cds_len = (trans_cds[-1][-1] - trans_cds[0][0]) + 1
        else:
            cds_len = 0
        cds_lengths.append(cds_len)

        # Exon lengths
        exons_lengths, introns_lengths = [[] for _ in range(2)]
        for (exon_1, exon_2) in trans_exons:
            exon_len = (exon_2 - exon_1) + 1
            exons_lengths.append(exon_len)

        # Track exon/intron lengths per transcript to generate table further downstream
        trans_exons_dt[trans_id].update(exons_lengths)

        # Intron lengths
        if trans_introns:
            for (intron_1, intron_2) in trans_introns:
                intron_len = (intron_2 - intron_1) + 1
                introns_lengths.append(intron_len)

            trans_introns_dt[trans_id].update(introns_lengths)

        # Memory efficient way to store large distributions for histograms
        exons_len_count = get_distribution_count(exons_lengths, exons_len_count, binwidth=1)
        if trans_introns:
            introns_len_count = get_distribution_count(introns_lengths, introns_len_count, binwidth=1)

    trans_len_count = get_distribution_count(trans_lengths, trans_len_count, binwidth=1)
    monoexon_len_count = get_distribution_count(monoexon_lengths, monoexon_len_count, binwidth=1)
    cds_len_count = get_distribution_count(cds_lengths, cds_len_count, binwidth=1)

    # Get quantiles from the lengths count
    lenght_data = [("Transcripts", trans_len_count), ("CDS", cds_len_count), ("Monoexonic Genes", monoexon_len_count),
                   ("Exons", exons_len_count), ("Introns", introns_len_count)]

    # Threshold quantiles must be present in quantile list
    lower_th, high_th = 1, 99
    perc = [0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99, 100]

    perc_data_dt, perc_ext_quants = [{} for _ in range(2)]
    for (cat_name, count_dt) in lenght_data:
        cat_perc = percentiles_from_counts(count_dt, percentiles_range=perc)
        cat_quants, cat_vals = list(zip(*cat_perc))

        perc_data_dt[cat_name] = cat_vals

        # Store lower and upper quantile values to generate tables with extreme cases
        low_ix = cat_quants.index(lower_th)
        lower_val = cat_vals[low_ix]

        high_ix = cat_quants.index(high_th)
        high_val = cat_vals[high_ix]

        perc_ext_quants[cat_name] = (lower_val, high_val)

    header = "Percentile\t" + "\t".join([str(v) for v in perc]) + "\n"
    lenghts_rows = [header]

    for cat_name, cat_vals in perc_data_dt.items():
        row = f"{cat_name}\t" + "\t".join([str(v) for v in cat_vals]) + "\n"
        lenghts_rows.append(row)

    # Write table
    simple_write_table(lenghts_rows, f"{outfile}_lengths_percentiles.tsv")

    # Get transcripts IDs in the lower and upper quantiles

    # Get the lower/upper values for each category
    low_quant_len, high_quant_len = perc_ext_quants["Transcripts"]
    low_quant_mono, high_quant_mono = perc_ext_quants["Monoexonic Genes"]
    low_quant_exon, high_quant_exon = perc_ext_quants["Exons"]
    low_quant_intron, high_quant_intron = perc_ext_quants["Introns"]

    outlier_len, outlier_mono, outlier_exon, outlier_intron = [[] for _ in range(4)]
    for t_id in gtf_obj.trans_exons_dt.keys():

        # Get the transcripts lengths values
        t_len = trans_lenght_dt[t_id]
        t_exon_lengths = sorted(trans_exons_dt[t_id])
        t_intron_lengths = sorted(trans_introns_dt[t_id])

        t_min_exon, t_max_exon = t_exon_lengths[0], t_exon_lengths[-1]
        if t_intron_lengths:
            t_min_intron, t_max_intron = t_intron_lengths[0], t_intron_lengths[-1]

        # Check for each category if the transcript is present below/above the outlier quantiles
        if t_len <= low_quant_len or t_len >= high_quant_len:
            t_dt = (t_id, t_len)
            outlier_len.append(t_dt)

        # Monoexon gene length
        if t_id in monoexon_trans:
            t_dt = (t_id, t_len)
            if t_len <= low_quant_mono or t_len >= high_quant_mono:
                outlier_mono.append(t_dt)

        # Exon length check
        if t_min_exon <= low_quant_exon:
            t_dt = (t_id, t_min_exon)
            outlier_exon.append(t_dt)
        elif t_max_exon >= high_quant_exon:
            t_dt = (t_id, t_max_exon)
            outlier_exon.append(t_dt)
        else:
            pass

        # Intron length check
        if t_intron_lengths:
            if t_min_intron <= low_quant_intron:
                t_dt = (t_id, t_min_intron)
                outlier_intron.append(t_dt)
            elif t_max_intron >= high_quant_intron:
                t_dt = (t_id, t_max_intron)
                outlier_intron.append(t_dt)
            else:
                pass

    length_id_data = [("Transcripts", outlier_len), ("Monoexonic Genes", outlier_mono), ("Exons", outlier_exon), ("Introns", outlier_intron)]

    for (cat_name, cat_vals) in length_id_data:
        cat_vals = sorted(cat_vals, key=lambda t_dt: t_dt[1])

        cat_rows = ["Transcript_ID\tFeature_length\n"]
        for (cat_id, len_val) in cat_vals:
            row = f"{cat_id}\t{len_val}\n"
            cat_rows.append(row)

        cat_tag = cat_name.lower().replace(" ", "_").replace("-", "_")
        simple_write_table(cat_rows, f"{outfile}_lengths_outliers_{cat_tag}.tsv")

    # Transcripts intron tables (potentially useful for analysis with compare with SJ reads table)
    trans_sj_rows = ["Scaffold\tTranscript_ID\tSJ_coordinates\n"]
    for t_id, t_introns in gtf_obj.trans_introns_dt.items():
        if t_introns:
            t_chr = gtf_obj.trans_chrom_dt[t_id]
            t_sj_lst = [f"{sj_1}-{sj_2}" for (sj_1, sj_2) in t_introns]

            row = f"{t_chr}\t{t_id}\t" + ",".join(t_sj_lst) + "\n"
            trans_sj_rows.append(row)

    simple_write_table(sorted(trans_sj_rows), f"{outfile}_IDs_splice_junction.tsv")

    return outpath

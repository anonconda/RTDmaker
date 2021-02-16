import os
import sys
import time
import glob
import shutil

from lib.tools.logger import logger
from lib.report.analyze_sj_data import generate_sj_report
from lib.report.analyze_annotation import generate_annotation_report
from lib.report.report_tools import simple_write_table, global_track_dt

from lib.filters.redundancy_QC import remove_redundant
from lib.filters.structural_QC import remove_chimeric_and_fragments
from lib.filters.GenLoc_QC import remove_ambiguous_location_models

from lib.filters.SJ_QC import identify_supported_sj, remove_unsupported_SJ_models
from lib.filters.abundance_QC import remove_low_abundant_via_quantification

from lib.tools.pad_annotation import pad_annotation
from lib.tools.reannotate_ids import generate_reannotated_gtf, generate_mapping_table
from lib.tools.quantify_transcripts import quantify_transcripts, generate_annotation_fasta
from lib.parsing.dir_parsing_tools import parse_annotation_table, generate_annotation_table


def is_empty(gtf_file):

    fl_size = os.stat(gtf_file).st_size
    if not fl_size:
        # Using "print + exit" instead of "sys.exit" so the error is printed into the logfile
        sys.exit(f'ERROR: Annotation file "{gtf_file}" is empty.')


def file_exist(outfile, skip_message=True):

    if os.path.exists(outfile):
        print(time.asctime(), f'File {outfile} already exist')
        if skip_message:
            print(time.asctime(), f'Skipping re-analysis and continuing with current existing file\n')
        else:
            print(time.asctime(), f'Overwriting file\n')
        return True
    else:
        return False


def remove_redundant_intermediary_files(dir_path, exceptions):

    intermediary_redundants = set()

    # Important: Tags "_temporary_union.gtf" and "_temporary_part" must be change if modified in the fuse_gtf()
    # and split_redundancy_removal() functions
    for temp_gtf in glob.glob(os.path.join(dir_path, f"**/*.gtf"), recursive=True):
        if temp_gtf.endswith("_temporary_union.gtf"):
            intermediary_redundants.add(temp_gtf)
        if "_temporary_part" in temp_gtf and temp_gtf not in exceptions:
            intermediary_redundants.add(temp_gtf)

    for temp_gtf in intermediary_redundants:
        os.remove(temp_gtf)


def main(ass_dir, ref_dir, sj_dir, sjreads_th, genome_fl, fastq_dir, abund_th,
         fragment_len, antisense_len, skip_set, add_set, size_th, keep_set, paths_dt, prefix, outname, logfile):

    # Arguments order and grouping:
    # srQC.main(args.ass_dir, args.ref_dir,                       # Annotations to analyze
    #           args.sj_dir, args.sjreads_th,                     # SJ-QC related args
    #           args.genome_fl, args.fastq_dir, args.abund_th,    # Abundance-QC related args
    #           args.len_th, args.antlen_th,                      # Structure-QC related args
    #           args.skip_set, args.add_set,                      # Category-selection related args
    #           args.size_th, args.keep_set,                      # Memory-related arguments
    #           paths_dt, args.prefix, args.outname, logfile)     # Output related args

    # Start logger
    logger(logfile)

    # Whether to print additional logging information or not
    verb = False

    # Track intermediary files to optionally remove them at the end
    intermediary_files = set()

    # Select whether write down removed files and/or keep intermediary files
    # Args if keep_set: 'intermediary', 'removed'

    # Parse the directories into tables (the table serve as a log of which files are process, tagged, and paired)
    assembly_table = generate_annotation_table(ass_dir, paths_dt["report"], "assemblies")
    references_table = generate_annotation_table(ref_dir, paths_dt["report"], "references")

    # Parse annotation tables into a list of files
    if assembly_table:
        gtf_files = parse_annotation_table(assembly_table)
    elif references_table:
        gtf_files = parse_annotation_table(references_table)
    else:
        sys.exit(f'ERROR: At least one table must be specified for the analysis')
    print("\n", end="")

    ############################################
    ## FIRST PART: ASSEMBLIES QUALITY CONTROL ##
    ############################################

    # Run this section if the assemblies are provided
    if ass_dir:
        # 1) Remove redundant transcripts
        # Redundancy removal already checks if file exist and avoid re-doing the analysis, no need to use file_exist()
        nonredundant_gtf = remove_redundant(gtf_files, paths_dt, f"{outname}_assemblies", size_th, to_keep=keep_set)
        is_empty(nonredundant_gtf)
        intermediary_files.add(nonredundant_gtf)

        # Removal of intermediary files (contain "_temporary_part" in their name or end with "_union.gtf") MUST be here
        if 'intermediary' not in keep_set:
            print(time.asctime(), f"Removing intermediary files")
            source_files = {nonredundant_gtf} | set(gtf_files)
            remove_redundant_intermediary_files(paths_dt['inter'], exceptions=source_files)
        print("\n", end="")

        # 2) Optional, Remove transcripts with uncertain strand information (e: unstranded, small antisense fragments)
        valid_loc_gtf = os.path.join(paths_dt['inter'], f"{outname}_assemblies_valid_loc.gtf")
        if not file_exist(valid_loc_gtf):
            valid_loc_gtf = remove_ambiguous_location_models(nonredundant_gtf, paths_dt, f"{outname}_assemblies",
                                                             to_add=add_set, to_keep=keep_set,
                                                             antisense_len=antisense_len, genome_fa=genome_fl)

            is_empty(valid_loc_gtf)
            intermediary_files.add(valid_loc_gtf)
            print("\n", end="")

        # 3) Remove models with poorly supported Splice Junctions
        supported_SJ_gtf = os.path.join(paths_dt['inter'], f"{outname}_supported_SJ.gtf")
        if "SJ" in skip_set:
            if file_exist(supported_SJ_gtf):
                pass
            else:
                print(time.asctime(), "WARNING: Skipping SJ-support QC analysis is not recommended.")
                supported_SJ_gtf = valid_loc_gtf
        else:
            # Check if the data to perform the SJ QC is available
            if not assembly_table or not sj_dir:
                sys.exit(f"ERROR: Data to perform the SJ Quality-Control is missing.")
            else:
                if file_exist(supported_SJ_gtf, skip_message=False):
                    pass

                SJ_sufolder = os.path.join(paths_dt["report"], "SJ_tables")
                sj_report = generate_sj_report(sj_dir, SJ_sufolder, outname)

                supported_sj_dt = identify_supported_sj(sj_dir, sjreads_th, paths_dt, outname)

                supported_SJ_gtf = remove_unsupported_SJ_models(valid_loc_gtf, supported_sj_dt, paths_dt, outname,
                                                                to_keep=keep_set)

                is_empty(supported_SJ_gtf)
                intermediary_files.add(supported_SJ_gtf)
                print("\n", end="")

        # 4) Remove transcripts poorly supported by reads (abundance)
        abundant_gtf = os.path.join(paths_dt["inter"], f"{outname}_abundant.gtf")
        if "quant" in skip_set:
            if file_exist(abundant_gtf):
                pass
            else:
                print(time.asctime(), "WARNING: Skipping abundance QC analysis is not recommended.")
                abundant_gtf = supported_SJ_gtf

        else:
            if not genome_fl or not fastq_dir:
                sys.exit(f"ERROR: Data to perform the transcripts quantification is missing.")
            else:
                if file_exist(abundant_gtf):
                    pass

                quant_dir = quantify_transcripts(supported_SJ_gtf, genome_fl, fastq_dir, paths_dt)
                print("\n", end="")

                abundant_gtf = remove_low_abundant_via_quantification(supported_SJ_gtf, quant_dir, paths_dt, outname,
                                                                      abund_th, to_keep=keep_set)

                is_empty(abundant_gtf)
                intermediary_files.add(abundant_gtf)
                print("\n", end="")

        accepted_gtf = abundant_gtf
        is_empty(accepted_gtf)
        intermediary_files.add(accepted_gtf)

        merge_files = [(accepted_gtf, f"_{outname}")]

    # If the Assembly QC is not executed, we need to initialize this variable
    else:
        accepted_gtf = None
        merge_files = []

    #############################################################
    ## SECOND PART: MERGE ASSEMBLIES WITH AVAILABLE REFERENCES ##
    #############################################################

    # The Merge QC is always executed after the Assembly QC, thus only check needed is to add reference info
    # Add references to the analysis (if provided)
    if references_table:
        merge_files += parse_annotation_table(references_table)

        print(time.asctime(), f"Integrating information from the references")
        if verb:
            for (gtf_fl, gtf_tag) in merge_files:
                print(time.asctime(), f"{gtf_fl}")
            print("\n", end="")
            # Sleep 1 second so that the lines are not redundant in logfile (redundant lines are deleted at the end)
            time.sleep(1)

    # 5) If there are multiple files to merge, we must check again for redundant & certain-strand transcripts

    # Previous version skipped the Redundancy and Ambiguous Location QC if only one reference was provided,
    # For current version we decided that these modules still need to performed even if only one reference is given
    # if len(merge_files) > 1:

    whole_anal = True
    if whole_anal:
        # IMPORTANT: the outfile generate by these methods ("remove_redundant", etc) MUST HAVE a different outname from
        # the files previously generate from the assembly_table, otherwise it will introduce a hidden bug!
        # The variable name also must be different, thus the 'ref_' prefix

        # Redundancy removal already checks if the file exist and avoid re-doing the analysis, no need to use file_exist
        ref_nonredundant_gtf = remove_redundant(merge_files, paths_dt, f"{outname}_references", size_th, to_keep=keep_set)

        is_empty(ref_nonredundant_gtf)
        intermediary_files.add(ref_nonredundant_gtf)
        print("\n", end="")

        # Removal of intermediary files (contain "_temporary_part" in their name or end with "_union.gtf") MUST be here
        if 'intermediary' not in keep_set:
            print(time.asctime(), f"Removing intermediary files")
            source_files = {ref_nonredundant_gtf} | set(merge_files) | set(gtf_files)
            remove_redundant_intermediary_files(paths_dt['inter'], exceptions=source_files)

        ref_valid_loc_gtf = os.path.join(paths_dt['inter'], f"{outname}_references_valid_loc.gtf")
        if not file_exist(ref_valid_loc_gtf):
            ref_valid_loc_gtf = remove_ambiguous_location_models(ref_nonredundant_gtf, paths_dt, f"{outname}_references",
                                                                 to_add=add_set, to_keep=keep_set,
                                                                 antisense_len=antisense_len, genome_fa=None)
            print("\n", end="")
    else:
        ref_valid_loc_gtf = merge_files[0][0]

    is_empty(ref_valid_loc_gtf)
    intermediary_files.add(ref_valid_loc_gtf)

    # 6) Check for overlaps (potentially chimeric models) and fragments
    refined_gtf = os.path.join(paths_dt["inter"], f"{outname}_to_reannotate.gtf")
    if not file_exist(refined_gtf):
        refined_gtf = remove_chimeric_and_fragments(ref_valid_loc_gtf, fragment_len, paths_dt, outname, logfile,
                                                         to_add=add_set, to_keep=keep_set)
        is_empty(refined_gtf)
        intermediary_files.add(refined_gtf)
        print("\n", end="")

    rtd_gtf = os.path.join(paths_dt["outpath"], f"{outname}.gtf")
    if not file_exist(rtd_gtf):
        rtd_gtf = generate_reannotated_gtf(refined_gtf, prefix, paths_dt, outname, logfile)
        is_empty(rtd_gtf)

    print(time.asctime(), f"RTD analysis completed!")
    print(time.asctime(), f"RTD location: {rtd_gtf}", flush=True)

    # Generate FASTA of RTD
    programs = {"gffread": shutil.which("gffread")}      # Absolute path og GFFREAD program
    generate_annotation_fasta(rtd_gtf, genome_fl, paths_dt["outpath"], programs)

    print("\n", end="")

    # Finish removing unused folders and unnecessary files
    if 'intermediary' not in keep_set:
        print(time.asctime(), f"Removing intermediary files")
        for inter_gtf in intermediary_files:
            os.remove(inter_gtf)

        # This remove the directory and all subfolders and files, INCLUDING the transcripts quantification analysis
        shutil.rmtree(paths_dt["inter"])

    if 'removed' not in keep_set:
        shutil.rmtree(paths_dt["removed"])


    #################################
    ## THIRD PART: GENERATE REPORT ##
    #################################

    print(time.asctime(), f"Starting RTD post-processing", flush=True)

    print(time.asctime(), f"Padding RTD file", flush=True)
    padded_gtf = os.path.join(paths_dt["outpath"], f"{outname}_padded.gtf")
    if not file_exist(padded_gtf):
        padded_gtf = pad_annotation(rtd_gtf, paths_dt, outname, logfile)
        is_empty(padded_gtf)

    # Generate FASTA of Padded version of RTD
    generate_annotation_fasta(padded_gtf, genome_fl, paths_dt["outpath"], programs)
    print("\n", end="")

    # Generate report tables of RTD
    rtd_sufolder = os.path.join(paths_dt["report"], "RTD_tables")
    generate_annotation_report(rtd_gtf, rtd_sufolder, outname, logfile)

    # Generate report tables of the references
    if references_table:
        references = parse_annotation_table(references_table)

        ref_sufolder = os.path.join(paths_dt["report"], "References_tables")
        for (ref_gtf, _) in references:
            ref_name = os.path.basename(ref_gtf).replace(".gtf", "")
            generate_annotation_report(ref_gtf, ref_sufolder, ref_name, logfile)

    # Generate report of the pipeline
    global_track_rows = ["Category\tFile\tNº of transcripts\n"]
    for fl, trans_n in global_track_dt.items():
        fl_cond, fl_name = fl.split("#")
        row = f"{fl_cond}\t{fl_name}\t{trans_n}\n"
        global_track_rows.append(row)

    global_track_file = os.path.join(paths_dt['report'], f"{outname}_pipeline_QC_steps_numbers.tsv")
    simple_write_table(global_track_rows, f"{global_track_file}")

    mapping_table = generate_mapping_table(rtd_gtf, paths_dt, outname)
    print("\n", end="")

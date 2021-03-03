import os
import sys
import time
import linecache
from collections import defaultdict
from itertools import islice, permutations

from lib.tools.logger import sizeof_fmt
from lib.tools.other_tools import get_transcript_length
from lib.parsing.gtf_object_tools import create_gtf_object, add_tag, write_gtf
from lib.tools.other_tools import sort_by_length, group_transcripts_by_overlap, get_start_end
from lib.report.report_tools import global_track_dt


def divide_by_size(files, size_th):

    if size_th < 1:
        sys.exit(f'The maximum file specified cannot be less than 1 Gb. Aborting.')

    # Convert size threshold from Gigabytes into bytes
    size_th *= 1000000000

    chunk, cum = [], 0
    for (fl, fl_tag) in files:
        fl_data = (fl, fl_tag)

        # File size is given in bytes
        fl_size = os.stat(fl).st_size
        cum += fl_size

        if cum <= size_th:
            chunk.append(fl_data)
        else:
            yield chunk
            chunk = [fl_data]
            cum = fl_size

    # To return the last element of file in list
    if chunk:
        yield chunk


def process_batch(gtf_files, paths_dt, outname, to_keep):

    verb = False
    if verb:
        print(time.asctime(), f"Batch files:")
        for (gtf_file, gtf_tag) in gtf_files:
            print(time.asctime(), f"{gtf_file}")
        print("\n")

    nonredundant_gtf = remove_redundant(gtf_files, paths_dt, outname, to_keep=to_keep)

    return nonredundant_gtf


def split_redundancy_removal(gtf_files, paths_dt, outname, size_th, to_keep):

    # If number of files to be processed is only 1, it wouldn't enter the "while" loop, thus this check to process it
    if len(gtf_files) == 1:
        outfile = process_batch(gtf_files, paths_dt, outname, to_keep)
        return outfile

    i = 1
    generated_files = []
    files_batches = gtf_files
    while len(files_batches) > 1:
        files_batches = list(divide_by_size(files_batches, size_th))

        n, temp = 1, []
        for batch in files_batches:
            outname = outname.split("_temporary_part")[0]
            outname += f"{i}"

            # If file already exist, ignore current iteration
            # Beware! the "predicted outfile" of iteration is the same name given by the "remove_redundants" method
            predicted_out = os.path.join(paths_dt["outpath"], f"{outname}_non_redundant.gtf")
            if os.path.exists(predicted_out):
                print(time.asctime(), f'File "{predicted_out}" already exist. Ignoring current iteration.')

                # Save current iteration output data for next iteration
                generated_files.append(predicted_out)
                file_data = (predicted_out, "")
                temp.append(file_data)
                n += 1
                i += 1
                continue

            print(time.asctime(), f"Processing batch {n} of {len(files_batches)}", flush=True)
            # Meld files and remove redundant transcripts
            batch_outfile = process_batch(batch, paths_dt, outname, to_keep)

            # Track files to remove them later
            generated_files.append(batch_outfile)

            # The tag is empty to avoid its re-annotation
            file_data = (batch_outfile, "")
            temp.append(file_data)
            n += 1
            i += 1

        # Perform the next iteration on the newly generated files
        files_batches = temp

    # Important! Sort files by creation time so that the file to keep is the last one!
    generated_files = sorted(generated_files, key=lambda f: os.path.getmtime(f))
    to_remove = generated_files[:-1]
    for fl in to_remove:
        os.remove(fl)

    outfile = generated_files[-1]

    return outfile


def fuse_gtf(gtf_files, outpath, outname):

    # If there is only one GTF file, there is no need to create an union file
    if len(gtf_files) == 1:
        # gtf_files format: [(gtf_file_path, gtf_tag_str)]
        return gtf_files[0][0]

    outfile = os.path.join(outpath, f"{outname}_temporary_union.gtf")

    print(time.asctime(), f'Merging files into a single annotation', flush=True)

    # If file already exist, overwrite it
    if os.path.exists(outfile):
        print(time.asctime(), f'File {outfile} already exist ({sizeof_fmt(os.stat(outfile).st_size)}).')

        overwrite = False
        if overwrite:
            print(time.asctime(), f'Overwriting file\n')
            with open(outfile, "w") as fh:
                pass
        else:
            print(time.asctime(), f'Keeping current file\n')
            return outfile

        # print(time.asctime(), f'File size: {sizeof_fmt(os.stat(outfile).st_size)}')

    tot = len(gtf_files)
    for n, (gtf_file, gtf_tag) in enumerate(gtf_files, 1):
        print(time.asctime(), f'Integrating file: {gtf_file} (file {n} out of {tot})')
        gtf_obj = create_gtf_object(gtf_file, gtf_tag)
        # Not using the "write_gtf" method because:
        # 1) the sorting (present in write_gtf) would unnecessarily slow the analysis
        # 2) this code contains the "add_tag" method necessary to make each line unique
        # 3) the "add_tag" method is unnecessary for write_gtf
        with open(outfile, "a+") as fh:
            for trans, lines_ix_list in gtf_obj.trans_gtf_lines_index.items():
                for ix in lines_ix_list:
                    line = linecache.getline(gtf_obj.gtf_path, ix)
                    # Allow to avoid re-annotation of tags
                    line = add_tag(line, gtf_tag)
                    fh.write(line)

        # print(time.asctime(), f'File updated ({outfile}), size: {sizeof_fmt(os.stat(outfile).st_size)})')

    return outfile


def identify_redundant_monoexons(gtf_obj, paths_dt):

    print(time.asctime(), f"Identifying redundant mono-exonic transcripts")

    # Track the redundant monoexons in a table, this file could be useful to choose which antisense to remove
    monoexon_redundant_dt = defaultdict(set)

    # Important! This method assumes that gtb_obj contains only mono-exonic transcripts!
    nonredundant_monoexons, redundant_monoexons = [set() for _ in range(2)]
    for chrom, strand_transcripts in gtf_obj.chrom_trans_dt.items():

        # Group transcripts by their overlap
        overlapping_transcripts = group_transcripts_by_overlap(gtf_obj, strand_transcripts)
        for overlap_group in overlapping_transcripts:
            # No need to check if group contains only mono-exons, as it is assumed by the method
            longest_trans = sort_by_length(gtf_obj, overlap_group)[0]
            redundants = [trans for trans in overlap_group if trans != longest_trans]

            nonredundant_monoexons.add(longest_trans)
            redundant_monoexons.update(redundants)

            monoexon_redundant_dt[longest_trans].update(redundants)

    return nonredundant_monoexons, redundant_monoexons


def identify_redundant_transcripts(gtf_obj):

    print(time.asctime(), "Identifying redundant transcripts with identical intron-coordinates")

    introns_trans_dt = defaultdict(set)
    for t_id, trans_introns in gtf_obj.trans_introns_dt.items():
        introns_trans_dt[tuple(trans_introns)].add(t_id)

    nonredundant_trans, redundant_trans, monoexonic_trans = [set() for _ in range(3)]
    for introns, trans_list in introns_trans_dt.items():
        # Mono-exonic transcripts must be processed separately
        if not introns:
            monoexonic_trans.update(trans_list)
            continue

        selected_trans = sort_by_length(gtf_obj, trans_list)[0]
        nonredundant_trans.add(selected_trans)

        redundants = [t_id for t_id in trans_list if t_id != selected_trans]
        redundant_trans.update(redundants)

    return nonredundant_trans, redundant_trans, monoexonic_trans


def slice_list(seq, n=2):

    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def is_subset(test_coordinates, query_coordinates):

    # Initial assumption
    flag = False

    # Slice list returns a tuple, thus we must convert this to make comparison
    # Convert here "test_coordinates" to make it more time efficient
    test_coordinates = tuple(test_coordinates)

    if len(test_coordinates) > len(query_coordinates):
        return flag

    for window in slice_list(query_coordinates, n=len(test_coordinates)):
        if window == test_coordinates:
            return True

    return flag


def identify_subset_redundant(gtf_obj, to_keep=None):

    print(time.asctime(), "Identifying redundant transcripts which coordinates are a subset of a larger model")

    redundant_subset = set()
    for chrom, strand_transcripts in gtf_obj.chrom_trans_dt.items():

        # Avoid analyzing transcript already known as redundant (to speed up analysis)
        if to_keep:
            strand_transcripts = [t_id for t_id in strand_transcripts if t_id in to_keep]

        # Group transcripts by their overlap
        overlapping_transcripts = group_transcripts_by_overlap(gtf_obj, strand_transcripts)
        for overlap_group in overlapping_transcripts:
            # Probably not necessary
            sorted_group = sorted(overlap_group, key=lambda t_id: get_transcript_length(gtf_obj, t_id), reverse=True)

            for trans_A, trans_B in permutations(sorted_group, 2):
                tA_introns = gtf_obj.trans_introns_dt[trans_A]
                tB_introns = gtf_obj.trans_introns_dt[trans_B]

                if tA_introns == tB_introns:
                    # For the moment ignore these transcripts, they are taking care of by the other function
                    continue

                if not tA_introns or not tB_introns:
                    # Ignore mono-exonics
                    continue

                if is_subset(tA_introns, tB_introns):
                    tA_start, tA_end = get_start_end(gtf_obj.trans_exons_dt[trans_A])
                    tB_start, tB_end = get_start_end(gtf_obj.trans_exons_dt[trans_B])

                    # Get the left/right coordinates from the intron to identify the exon PAIRS in trans_B
                    ex_left = tA_introns[0][0] - 1
                    ex_right = tA_introns[-1][-1] + 1

                    boundary_exon_left, boundary_exon_right = (), ()
                    for exon in gtf_obj.trans_exons_dt[trans_B]:
                        if ex_left in exon:
                            boundary_exon_left = exon
                        if ex_right in exon:
                            boundary_exon_right = exon

                    if not boundary_exon_left or not boundary_exon_right:
                        continue

                    left_flag, right_flag = False, False
                    if boundary_exon_left[0] <= tA_start <= boundary_exon_left[1]:
                        left_flag = True

                    if boundary_exon_right[0] <= tA_end <= boundary_exon_right[1]:
                        right_flag = True

                    if left_flag and right_flag:
                        redundant_subset.add(trans_A)

    return redundant_subset


def remove_redundant(gtf_files, paths_dt, outname, size_th=None, verb=True, to_keep=None, disable=False):

    # Important! The "split_analysis" predict this outfile to avoid running pre-calculated analysis
    # If you change the format here, you must also change it on the "split_analysis" method

    if not to_keep:
        to_keep = set()

    # This function requires the 'intermediary' files to do its analysis, this are removed later on srQC_main
    to_keep.add('intermediary')

    accepted_outname = f"{outname}_non_redundant.gtf"

    predicted_out = os.path.join(paths_dt['inter'], accepted_outname)

    # Allow to ignore processing the files that already exist (to speed up analysis of already processed large datasets)
    if os.path.exists(predicted_out):
        print(time.asctime(), f'File "{predicted_out}" already exist ({sizeof_fmt(os.stat(predicted_out).st_size)})')

        # If file already exist, overwrite it
        overwrite = False
        if overwrite:
            print(time.asctime(), f'Overwriting file')
            with open(predicted_out, "w") as fh:
                pass
        else:
            print(time.asctime(), f'Keeping current file')
            return predicted_out

    if verb:
        print(time.asctime(), f"Removing redundant transcripts")
        for (gtf_fl, gtf_tag) in gtf_files:
            # Print names into logfile (disabled)
            pass

    if size_th:
        nonredundant_gtf = split_redundancy_removal(gtf_files, paths_dt, outname, size_th, to_keep)
        return nonredundant_gtf

    # Put all transcripts models into a single file, this function also make the IDs unique for each file
    union_gtf = fuse_gtf(gtf_files, paths_dt["inter"], outname)

    gtf_obj = create_gtf_object(union_gtf)

    # Identify multi-exonic redundant
    nonredundant_trans, redundant_trans, monoexonic_trans = identify_redundant_transcripts(gtf_obj)

    redundant_by_subset = identify_subset_redundant(gtf_obj, to_keep=nonredundant_trans)

    # Remove the subset models from previous accepted "non-redundant"
    redundant_by_subset = redundant_by_subset - redundant_trans
    nonredundant_trans = nonredundant_trans - (redundant_trans | redundant_by_subset)

    # Redundancy is defined by intronic-structure, thus mono-exons must be processed separately
    gtf_obj_monoexons = create_gtf_object(union_gtf, to_keep=monoexonic_trans)

    nonredundant_monoexons, redundant_monoexons = identify_redundant_monoexons(gtf_obj_monoexons, paths_dt)

    # Important: Re-introduce the non-redundant mono-exons into the non-redundant output file!
    nonredundant_trans = nonredundant_trans | nonredundant_monoexons

    # Optionally keep or remove redundant transcripts by intron-coordinates, i.e. when analyzing data with PacBio data
    if disable:
        nonredundant_trans.update(redundant_trans | redundant_by_subset)
        redundant_trans, redundant_by_subset = set(), set()

    # Redundant transcript
    tot_redundants = redundant_trans | redundant_by_subset | redundant_monoexons

    outfile = write_gtf(gtf_obj, nonredundant_trans, paths_dt['inter'], accepted_outname)

    # Important: to keep SOURCE files, do NOT remove original files from inside function, use external function instead
    if 'removed' in to_keep:
        _ = write_gtf(gtf_obj, redundant_trans, paths_dt['removed'], f"{outname}_redundant_intron_coord_identical.gtf")
        _ = write_gtf(gtf_obj, redundant_by_subset, paths_dt['removed'], f"{outname}_redundant_intron_coord_subset.gtf")
        _ = write_gtf(gtf_obj, redundant_monoexons, paths_dt['removed'], f"{outname}_redundant_monoexons.gtf")

    # Track numbers of accepted / removed
    global_track_dt[f"Removed#{outname}_redundant_intron_coord_identical.gtf"] = len(redundant_trans)
    global_track_dt[f"Removed#{outname}_redundant_intron_coord_subset.gtf"] = len(redundant_by_subset)
    global_track_dt[f"Removed#{outname}_redundant_monoexons.gtf"] = len(redundant_monoexons)

    global_track_dt[f"Accepted#{accepted_outname}"] = len(nonredundant_trans)

    return outfile



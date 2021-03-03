import os
import sys
import time
import math
from lib.tools.other_tools import flat
from collections import OrderedDict, defaultdict

# This dictionary is to track the number of accepted and removed models at each Quality Control step
global_track_dt = OrderedDict()


def percentiles_from_counts(counts_dt, percentiles_range=None):
    """Returns [(percentile, value)] with nearest rank percentiles.
    Percentile 0: <min_value>, 100: <max_value>.
    counts_dt: { <value>: <count> }
    percentiles_range: iterable for percentiles to calculate; 0 <= ~ <= 100
    Source: https://stackoverflow.com/questions/25070086/percentiles-from-counts-of-values
    """

    if percentiles_range is None:
        percentiles_range = [0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99, 100]

    # This handle the input if the case under analysis is completely absent in the annotation (count_dt is empty)
    if not counts_dt:
        vals = [0 for _ in range(len(percentiles_range))]
        percentiles = list(zip(percentiles_range, vals))
        return percentiles

    assert all(0 <= p <= 100 for p in percentiles_range)

    percentiles = []
    num = sum(counts_dt.values())
    counts = sorted(counts_dt.items())
    curr_counts_pos = 0  # current position in counts
    curr_pos = counts[0][1]  # sum of frequencies up to curr_counts_pos
    for p in sorted(percentiles_range):
        if p < 100:
            percentile_pos = p / 100.0 * num
            while curr_pos <= percentile_pos and curr_counts_pos < len(counts):
                curr_counts_pos += 1
                curr_pos += counts[curr_counts_pos][1]
            percentiles.append((p, counts[curr_counts_pos][0]))
        else:
            percentiles.append((p, counts[-1][0]))

    return percentiles


def percentile(N, percent, key=lambda x: x):

    # Ensure the user is using the correct input
    assert 0.0 <= percent <= 1.0

    if not N:
        return None
    # This sort was added by me to avoid possible errors if N is not sorted
    N = sorted(N)
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        try:
            return key(N[int(k)])
        except IndexError as e:
            sys.exit(f"Error while calculating percentile. Please check logfile for full error traceback.\n{e}")

    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1


def get_per(n, tot):

    try:
        return round((n/tot)*100, 1)
    except ZeroDivisionError:
        return 0.0


def get_location(f_id, gtf_obj, feature="trans"):

    if feature == "trans":
        try:
            chrom = gtf_obj.trans_chrom_dt[f_id][:-1]
            t_exons = flat(gtf_obj.trans_exons_dt[f_id])
            coord = (min(t_exons), max(t_exons))
        except KeyError:
            return "-"

    elif feature == "gene":
        try:
            t_id = sorted(gtf_obj.gene_trans_dt[f_id])[0]
            chrom = gtf_obj.trans_chrom_dt[t_id][:-1]
            coord = gtf_obj.gene_coords_dt[f_id]
        except KeyError:
            return "-"

    else:
        sys.exit(f'Feature type must be either "trans" or "gene", not "{feature}"')

    coord_txt = f'{chrom}:{coord[0]}-{coord[1]}'

    return coord_txt


def get_sort_key(t_id, gtf_obj):

    # Sort transcripts by the key (Chrom, Gene_ID, Trans_ID, Leftmost_coordinate)
    return (gtf_obj.trans_chrom_dt[t_id], gtf_obj.trans_gene_dt[t_id], t_id, gtf_obj.trans_exons_dt[t_id][0][0])


def get_distribution_count(values, counts_dt, binwidth=1):

    # Memory efficient way to store large distributions for histogram:
    # 1) data is rounded up to the nearest bindwith
    # 2) the rounded values is used as key of a dictionary tracking the count of how many times the value is seen
    # Source: https://stackoverflow.com/questions/2464871/numpy-histogram-of-large-arrays
    for val in values:
        binname = int(math.floor(val / binwidth) * binwidth)
        if binname not in counts_dt:
            counts_dt[binname] = 0
        counts_dt[binname] += 1

    return counts_dt


def write_table(transcripts, gtf_obj, outpath, outname):

    # Unused function TODO check if it is safe to delete

    # Create an output folder if it doesnt exist
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    outfile = os.path.join(outpath, outname)

    # Sort transcripts by the key (Chrom, Gene_ID, Trans_ID, Leftmost_coordinate)
    get_sort_key = lambda t_id: (gtf_obj.trans_chrom_dt[t_id], gtf_obj.trans_gene_dt[t_id], t_id,
                                 gtf_obj.trans_exons_dt[t_id][0][0])

    sorted_transcripts = sorted(transcripts, key=lambda t_id: get_sort_key(t_id))

    print(time.asctime(), f'Writing table: {outfile}')
    with open(outfile, "w+") as fh:
        # Header
        fh.write(f'Location,Transcript_ID\n')
        for trans_id in sorted_transcripts:
            if trans_id in transcripts:
                trans_loc = get_location(trans_id, gtf_obj, "trans")
                row = f'{trans_loc},{trans_id}\n'
                fh.write(row)


def simple_write_table(rows_list, outfile):

    with open(outfile, "w+") as fh:
        for row in rows_list:
            fh.write(row)

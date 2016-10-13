"""
A script to collect initialization task statuses counts and compute time stats
among the children or grandchildren of some containers.
"""

from __future__ import division
from genestack_client import make_connection_parser, get_connection, FilesUtil
from collections import Counter


def chunks(l, n):
    """Yield successive n-sized chunks from l (except from last chunk which may be smaller)."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

STATUSES = ['NOT_STARTED', 'IN_PROGRESS', 'FAILED', 'COMPLETE']


def compute_init_duration(time_map):
    if 'initializationStart' not in time_map or 'initializationEnd' not in time_map:
        return 0
    start = time_map['initializationStart']/1000
    end = time_map['initializationEnd']/1000
    return end - start


def get_task_counts_and_compute_time(files):
    counts = Counter()
    compute_time = 0
    computed_files = 0
    for batch in chunks(files, 200):
        infos = fu.get_infos(batch)
        for info in infos:
            status = info['initializationStatus']['id']
            counts[status] += 1
            t = compute_init_duration(info['time'])
            if t > 0:
                computed_files += 1
                compute_time += t
    return counts, compute_time, computed_files


if __name__ == "__main__":
    p = make_connection_parser()
    p.add_argument('folders', help='Container(s) accession(s)', nargs="+")
    p.add_argument('-rr', action="store_true", help='Flag to collect all grandchildren of the provided containers,'
                                 ' instead of children')
    p.add_argument('--name', action='store_true', help='Display folder name')
    args = p.parse_args()

    connection = get_connection(args)
    fu = FilesUtil(connection)

    all_folders = []
    if not args.rr:
        all_folders = args.folders
    else:
        for fd in args.folders:
            all_folders += fu.get_file_children(fd)

    metainfo_keys = ['accession', 'computed', 'total_time(h)', 'avg_time(h)']
    if args.name:
        metainfo_keys = ['name'] + metainfo_keys
    print "\t".join(metainfo_keys + STATUSES)

    for fd in all_folders:
        counts, time, computed = get_task_counts_and_compute_time(fu.get_file_children(fd))
        fd_info = fu.get_infos([fd])[0]
        custom_info = {
            'accession': fd_info['accession'],
            'name': fd_info['name'],
            'computed': computed,
            'total_time(h)': round(time/3600, 2),
            'avg_time(h)': round(time/(3600*computed), 2) if computed > 0 else "NA"
        }
        metainfo_values = map(str, [custom_info[k] for k in metainfo_keys])
        counts_values = map(str, [counts[key] for key in STATUSES])
        print "\t".join(metainfo_values + counts_values)

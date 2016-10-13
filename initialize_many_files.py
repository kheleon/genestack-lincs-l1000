"""
A script to monitor the initialization of a container's children.
The script will keep initializing files in the container while there are still files that
can be initialized (NOT_STARTED or FAILED).
You can also enforce to have a maximum number of tasks active (starting or running) using -m
"""

from genestack_client import make_connection_parser, get_connection, FilesUtil, FileInitializer, Application
from collections import Counter
import time
import sys

BATCH_SIZE = 10
WAIT_TIME = 60


def chunks(l, n):
    """Yield successive n-sized chunks from l (except from last chunk which may be smaller)."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


class TaskManager(Application):
    APPLICATION_ID = 'genestack/taskmanager'

    def _get_tasks(self, status):
        return self.invoke('getTasks',
                           {'status': status}, "file", "asc", 0, 0)['count']

    def get_active_tasks(self):
        return self._get_tasks("RUNNING") + self._get_tasks("STARTING")

p = make_connection_parser()
p.add_argument('folder', help='Accession of the container whose children will be initialized')
p.add_argument('-i', '--info', help='Show initialization status counts and quit', action='store_true')
p.add_argument('-m', '--max', help='Set a maximum number of simultaneously active tasks', type=int, default=-1)
args = p.parse_args()

print "Connecting to Genestack..."
connection = get_connection(args)
fu = FilesUtil(connection)
fi = FileInitializer(connection)
tm = TaskManager(connection)

print "Collecting container children..."
children = fu.collect_initializable_files_in_container(args.folder)
nb_children = len(children)


def print_counts(counter):
    print "Initialization statuses:"
    for k,v in counter.iteritems():
        print " * %s: %d" % (k,v)


def get_initializable_files(children):
    print "Collecting initializable files..."
    counts = Counter()
    initializable = []
    for batch in chunks(children, 200):
        infos = fu.get_infos(batch)
        for info in infos:
            status = info['initializationStatus']['id']
            if status == "NOT_STARTED" or status == "FAILED":
                initializable.append(info['accession'])
            counts[status] += 1
    return counts, initializable


counts, initializable = get_initializable_files(children)
print_counts(counts)

if args.info:
    sys.exit(0)

while counts['COMPLETE'] < sum(counts.values()):
    while len(initializable) == 0 and counts['COMPLETE'] < sum(counts.values()):
        print "No files can be initialized yet. Waiting..."
        time.sleep(60)
        counts, initializable = get_initializable_files(children)
        print_counts(counts)

    initialized_count = 0
    for batch in chunks(initializable, BATCH_SIZE):
        if args.max > -1:
            time.sleep(5)
            active = tm.get_active_tasks()
            while active >= args.max - BATCH_SIZE:
                print "\n%d active tasks. Waiting for the queue to clear up..." % active
                time.sleep(WAIT_TIME)
                active = tm.get_active_tasks()

        initialized_count += len(batch)
        sys.stdout.write("\rInitializing... (%d/%d)" % (initialized_count, len(initializable)))
        sys.stdout.flush()
        fi.initialize(batch)

    print ""
    counts, initializable = get_initializable_files(children)
    print_counts(counts)

print "All files were initialized successfully."


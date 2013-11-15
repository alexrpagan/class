import os
import argparse
import sys
import glob

from Bio.Blast.Applications import NcbiblastnCommandline
from multiprocessing.pool import ThreadPool
from utils import validdir

QUERYFMT = "*.txt"
CMD = "blastn"
OUTFMT = 7

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("querydir", help="path to directory of query files", type=str)
    parser.add_argument("dbdir",    help="path to directory of databases", type=str)
    parser.add_argument("outdir",   help="path to directory of output files", type=str)
    parser.add_argument("-e", "--evalue", help="evalue for BLAST", type=float, default=1e-1)
    parser.add_argument("-t", "--title", help="base title for outfile", type=str, default="out")
    args = parser.parse_args()

    blastpath = None
    if 'BLASTPATH' in os.environ:
        blastpath = os.environ['BLASTPATH']
    if not validdir(blastpath, "Invalid $BLASTPATH"):
        return 1
    if not validdir(args.dbdir, "Missing database directory"):
        return 1
    if not validdir(args.querydir, "Missing or invalid query file directory"):
        return 1

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    blastcmd = os.path.join(blastpath, CMD)
    queries = get_query_files(args.querydir)
    dbs = get_db_names(args.dbdir)

    jobs = []
    for i, query in enumerate(queries):
        for j, db in enumerate(dbs):
            outfilename = "%s_%d_%s" % (args.title, i, j)
            job = NcbiblastnCommandline(blastcmd
                , query=query
                , db=db
                , evalue=args.evalue
                , outfmt=OUTFMT
                , out=os.path.join(args.outdir, outfilename))
            jobs.append(job)

    pool = ThreadPool()
    pool.map(get_job_executor(jobs), range(len(jobs)))

    return 0

def get_job_executor(jobs):
    def do_job(i):
        job = jobs[i]
        outfile = job.out
        job()
        process_output_file(i, outfile)
    return do_job

def process_output_file(job_number, filepath):
    print "Processed output file for job %d" % job_number

def get_query_files(querydir):
    return glob.glob(os.path.join(querydir, QUERYFMT))

def get_db_names(dbdir):
    dbnames = set()
    files = glob.glob(os.path.join(dbdir, "*"))
    for f in files:
        base = os.path.basename(f)
        dirname = os.path.dirname(f)
        s = base.split('.')
        if len(s) == 2:
            if s[0] not in dbnames:
                dbnames.add(os.path.join(dirname, s[0]))
    return list(dbnames)

if __name__ == "__main__":
    sys.exit(main())

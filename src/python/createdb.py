import argparse
import os
import sys
import glob
import subprocess

from utils import chunks, validdir

INFORMAT = '*.fa'
DBTYPE = 'nucl'
BIN = 'makeblastdb'

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("indir", help="path to directory of input files", type=str)
    parser.add_argument("outdir", help="path to directory of database files", type=str)
    parser.add_argument("-t", "--title", help="the base of the title for the blastdb", type=str, default="blastdb")
    parser.add_argument("-p", "--partitions", help="number of files to partition database into", type=int, default=1)
    parser.add_argument("--minpar", help="use the partition size", type=bool)
    args = parser.parse_args()

    blastpath = None
    if 'BLASTPATH' in os.environ:
        blastpath = os.environ['BLASTPATH']
    if not validdir(blastpath, "Invalid $BLASTPATH"):
        return 1
    if not validdir(args.indir, "Input directory does not exist"):
        return 1

    infiles = glob.glob(os.path.join(args.indir, INFORMAT))
    if len(infiles) == 0:
        sys.stderr.write("No valid input files")
        return 1

    num_partitions = args.partitions
    if args.minpar:
        num_partitions = len(infiles)

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    processes = []
    partitions = list(chunks(infiles, len(infiles)/num_partitions))
    for i, partition in zip(range(len(partitions)), partitions):
        processes.append(createdb(blastpath, i, partition, args.outdir, args.title))

    for proc in processes:
        out, err = proc.communicate()
        print out

    return 0

def createdb(blastpath, run_num, infiles, outdir, title):
    infilestr  = ' '.join(infiles)
    titlestr   = "%s%d" % (title, run_num)
    outfilestr = os.path.join(outdir, titlestr)
    binpath    = os.path.join(blastpath, BIN)
    cmd = [binpath, '-in', infilestr, '-out', outfilestr, '-dbtype', DBTYPE, '-title', titlestr]
    return subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=False)

if __name__ == "__main__":
    sys.exit(main())


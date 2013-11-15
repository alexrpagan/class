import sys
import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="variant db file", type=str)
    args = parser.parse_args()
    tag = get_chr_tag(os.path.basename(args.filename))
    tagsize = len(tag)
    with open(args.filename, 'r') as f:
        for line in f:
            pos = line.split()[1]
            line = line.rstrip()
            fmt = "*4\r\n$4\r\nZADD\r\n$%d\r\n%s\r\n$%d\r\n%s\r\n$%d\r\n%s\r\n"
            sys.stdout.write(fmt % (tagsize, tag, len(pos), pos, len(line), line))
    return 0

def get_chr_tag(basename):
    tag = ''
    for part in basename.split('.'):
        if 'chr' in part:
            tag += part
        elif 'vcf' in part:
            tag += ('-' + part)
    return tag

if __name__ == "__main__":
    sys.exit(main())
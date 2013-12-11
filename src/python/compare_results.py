import sys
import collections

FUZZ = 1000
BENCHLOC='/datadisk1/apagan/benchresults/'

def main():
    if len(sys.argv) < 2:
        print "Too few arguments!"
        return

    size = sys.argv[1]

    class_ignores = set()
    # get queries to ignore
    with open(BENCHLOC + 'class_' + size + '/err') as f:
        for line in f:
            line = line.strip()
            if 'Skipping!' in line:
                class_ignores.add(int(line.split()[4]))

    cablast_ignores = set()
    with open(BENCHLOC + 'cablast_' + size + '/err') as f:
        for line in f:
            line = line.strip()
            if 'skipping' in line:
                cablast_ignores.add(int(line.split()[3][:-1]))

    # print class_ignores & cablast_ignores
    # print class_ignores - cablast_ignores
    # print cablast_ignores - class_ignores
    ignore = class_ignores | cablast_ignores

    blast_results = collections.defaultdict(list)
    with open(BENCHLOC + 'blast_' + size + '/out') as f:
        query_idx = -1
        for line in f:
            line = line.strip()
            if 'BLASTN 2.2.28+' in line:
                query_idx += 1
            else:
                if line[0] != "#":
                    if query_idx not in ignore:
                        blast_results[query_idx].append(line.split('\t'))

    class_results = collections.defaultdict(list)
    with open(BENCHLOC + 'class_' + size + '/out') as f:
        query_idx = -1
        for line in f:
            line = line.strip()
            if line[0] == '#':
                query_idx = int(line.split()[2])
            else:
                class_results[query_idx].append(line.split('\t'))


    cablast_results = collections.defaultdict(list)
    with open(BENCHLOC + 'cablast_' + size + '/out') as f:
        query_idx = -1
        for line in f:
            line = line.strip()
            if line[0] == '#':
                query_idx = int(line.split()[2])
            else:
                cablast_results[query_idx].append(line.split('\t'))

    hits = 0
    class_misses   = 0
    cablast_misses = 0

    for query in blast_results:
        for hit in blast_results[query]:
            hits += 1
            indiv = hit[1]
            start = int(hit[8])
            end   = int(hit[9])
            found = False
            for class_hit in class_results[query]:
                if convert_notation(class_hit[0]) == indiv:
                    class_start = int(class_hit[1])
                    class_end   = int(class_hit[2])
                    if (class_start - FUZZ <= end and class_end + FUZZ >= start):
                        found = True
                        break
            if not found:
                class_misses += 1
    class_recall = 1 - (float(class_misses) / hits)

    class_presc_miss = 0
    class_hits = 0
    for query in blast_results:
        for class_hit in class_results[query]:
            class_hits += 1
            indiv = convert_notation(class_hit[0])
            class_start = int(class_hit[1])
            class_end   = int(class_hit[2])
            found = False
            for hit in blast_results[query]:
                if hit[1] == indiv:
                    start = int(hit[8])
                    end   = int(hit[9])
                if (class_start - FUZZ <= end and class_end + FUZZ >= start):
                    found = True
                    break
            if not found:
                class_presc_miss += 1
    class_presc =  1 - (float(class_presc_miss) / class_hits)

    print "CLASS   presc: %f, recall: %f, fscore: %f " %  (class_presc, class_recall, fscore(class_presc, class_recall))

    for query in blast_results:
        for hit in blast_results[query]:
            indiv = hit[1]
            start = int(hit[8])
            end   = int(hit[9])
            found = False
            for cablast_hit in cablast_results[query]:
                if cablast_hit[0] == indiv:
                    cablast_start = int(cablast_hit[1])
                    cablast_end   = int(cablast_hit[2])
                    if (cablast_start - FUZZ <= end and cablast_end + FUZZ >= start):
                        found = True
                        break
            if not found:
                cablast_misses += 1
    cablast_recall =  1 - (float(cablast_misses) / hits)

    cablast_presc_miss = 0
    cablast_hits = 0
    for query in blast_results:
        for cablast_hit in cablast_results[query]:
            cablast_hits += 1
            indiv = cablast_hit[0]
            cablast_start = int(cablast_hit[1])
            cablast_end   = int(cablast_hit[2])
            found = False
            for hit in blast_results[query]:
                if hit[1] == indiv:
                    start = int(hit[8])
                    end   = int(hit[9])
                if (cablast_start - FUZZ <= end and cablast_end + FUZZ >= start):
                    found = True
                    break
            if not found:
                cablast_presc_miss += 1
    cablast_presc = 1 - (float(cablast_presc_miss) / cablast_hits)

    print "CaBLAST presc: %f, recall: %f, fscore: %f " %  (cablast_presc, cablast_recall, fscore(cablast_presc, cablast_recall))

def fscore(p,r):
    return (r * p) / (r + p)

def convert_notation(indiv):
    parts = indiv.split('-')
    replace = '1'
    if parts[1] == '1':
        replace = '0'
    return '-'.join([parts[0], replace])

if __name__ == "__main__":
    main()
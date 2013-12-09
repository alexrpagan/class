import random

REFERENCE   = '/datadisk1/apagan/outbench/ref/chr22.fa'
NUM_QUERIES = 5000
MIN_LENGTH  = 5
MAX_LENGTH  = 200
SEED        = 1

#random.seed(SEED)

BASES = ['A', 'C', 'G', 'T']

def main():
    queries = []
    origins = []
    ref = ""
    seen_defline = False
    with open(REFERENCE, 'r') as f:
        for line in f:
            if not seen_defline:
                seen_defline = True
                continue
            ref += line.strip()

    for _ in range(NUM_QUERIES):
        start   = random.randint(0, len(ref) - MAX_LENGTH - 1)
        length  = random.randint(MIN_LENGTH, MAX_LENGTH)
        replace = random.random()
        ref_seq = ref[start:start+length]
        query = ""
        for c in ref_seq:
            if random.random() < replace:
                query += random.choice(BASES)
            else:
                query += c
        queries.append(query)
        origins.append((start, replace,))

    for q, o in zip(queries, origins):
        print ">%s|origin: %d|replace:%f" % (REFERENCE, o[0], o[1])
        print q

if __name__ == "__main__":
    main()

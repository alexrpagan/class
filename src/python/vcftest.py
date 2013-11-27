import vcf

def main():
    i = 0
    vcf_reader = vcf.Reader(open('/home/apagan/code/tgc/chr22/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf'))
    for record in vcf_reader:
        i = 1 - i

if __name__ == "__main__":
    main()
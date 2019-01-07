class CHROMS_VCF:
    '''
    Parse VCF file and split records into different Chromosomes.
    '''
    def __init__(self, file_name):
        
        self.parse_vcf_file(file_name)

        
    def parse_vcf_file(self, file_name):
        
        self.CHROMS = {}
        with open(file_name) as vcf_handle:
            for vcf_record in vcf_handle:
                if vcf_record.startswith('##'):
                    continue
                vcf_field_content = vcf_record.split('\t')
                CHROM, POS, ID, REF, ALT, *OTHERS = vcf_field_content
                if not self.is_point_mutation(REF, ALT):
                    continue
                if not CHROM.isdigit():
                    continue
                CHROM = int(CHROM)
                if CHROM in self.CHROMS:
                    self.CHROMS[CHROM].append((int(POS), REF, ALT))
                else:
                    self.CHROMS[CHROM] = [(int(POS), REF, ALT)]
        for key in self.CHROMS:
            self.CHROMS[key].sort()
            
            
    def is_point_mutation(self, REF, ALT):
        
        return len(REF) == 1 and len(ALT) == 1
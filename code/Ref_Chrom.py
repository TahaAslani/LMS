import gzip
from Bio import SeqIO

class Ref_Chrom:
    '''
    Load Human Chromosome from reference fasta file by index.
    '''
    def __init__(self, index):
        
        seq_list = []
        with gzip.open("../ref/GrCH37_lite/chr" + str(index) +".fa.gz", "rt") as handle:
            for seq in SeqIO.parse(handle, "fasta"):
                seq_list.append(seq)
                
        assert len(seq_list) == 1, '''ERROR: The number of sequences should equal to 1 
                       (it is likely that the reference genomes are not assembled)'''
        
        self.seq = str(seq_list[0].seq)
        self.len = len(self.seq)
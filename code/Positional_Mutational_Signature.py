from Interval import Interval
from CHROMS_VCF import CHROMS_VCF
from Ref_Chrom import Ref_Chrom

import numpy as np
import matplotlib
from matplotlib import pyplot as plt

class Positional_Mutational_Signature:
    '''
    Computing positional mutational signature and overall signature
    Arguments:
           WINDOW_SIZE -- int, the legnth of window per positional ms
           subject_vcf_path -- vcf file path per subject
           filter_size -- compute pms when vcf records in the window is greater than this threshold
           mode -- 0: ignore A -> (T, G, C) and G -> (A, T, C) mutations
                   1: convert A and G type to C and G type mutations
    Result:
         pms -- a instance of the class
    '''
    def __init__(self, WINDOW_SIZE = 10**6, subject_file_path = '../Texas_vcfs/wheeljack1.vcf', filter_size = 96, mode = '1'):
        self.subject_file_path = subject_file_path
        self.mode = mode
        self.filter_size = filter_size
        self.bad_records = 0
        self.overall_mutational_signature = {'ACA-A':0,'ACC-A':0,'ACG-A':0,'ACT-A':0,'CCA-A':0,
                                             'CCC-A':0,'CCG-A':0,'CCT-A':0,'GCA-A':0,'GCC-A':0,
                                             'GCG-A':0,'GCT-A':0,'TCA-A':0,'TCC-A':0,'TCG-A':0,
                                             'TCT-A':0,'ACA-G':0,'ACC-G':0,'ACG-G':0,'ACT-G':0,
                                             'CCA-G':0,'CCC-G':0,'CCG-G':0,'CCT-G':0,'GCA-G':0,
                                             'GCC-G':0,'GCG-G':0,'GCT-G':0,'TCA-G':0,'TCC-G':0,
                                             'TCG-G':0,'TCT-G':0,'ACA-T':0,'ACC-T':0,'ACG-T':0,
                                             'ACT-T':0,'CCA-T':0,'CCC-T':0,'CCG-T':0,'CCT-T':0,
                                             'GCA-T':0,'GCC-T':0,'GCG-T':0,'GCT-T':0,'TCA-T':0,
                                             'TCC-T':0,'TCG-T':0,'TCT-T':0,'ATA-A':0,'ATC-A':0,
                                             'ATG-A':0,'ATT-A':0,'CTA-A':0,'CTC-A':0,'CTG-A':0,
                                             'CTT-A':0,'GTA-A':0,'GTC-A':0,'GTG-A':0,'GTT-A':0,
                                             'TTA-A':0,'TTC-A':0,'TTG-A':0,'TTT-A':0,'ATA-C':0,
                                             'ATC-C':0,'ATG-C':0,'ATT-C':0,'CTA-C':0,'CTC-C':0,
                                             'CTG-C':0,'CTT-C':0,'GTA-C':0,'GTC-C':0,'GTG-C':0,
                                             'GTT-C':0,'TTA-C':0,'TTC-C':0,'TTG-C':0,'TTT-C':0,
                                             'ATA-G':0,'ATC-G':0,'ATG-G':0,'ATT-G':0,'CTA-G':0,
                                             'CTC-G':0,'CTG-G':0,'CTT-G':0,'GTA-G':0,'GTC-G':0,
                                             'GTG-G':0,'GTT-G':0,'TTA-G':0,'TTC-G':0,'TTG-G':0,'TTT-G':0}

        self.pms_dict = {}
        self.pms_burden = {}

        vcf_subject = CHROMS_VCF(subject_file_path)
        
        for chroms_index in range(1, 23):
            print('[PROGRESS] processing Chromosome #%d' % chroms_index)
            # load human reference chromosome
            ref_chrom = Ref_Chrom(chroms_index)
            # generate window interal objects
            intervals = self.window_generator(ref_chrom.len, WINDOW_SIZE)
            # split vcf records into windows
            if chroms_index not in vcf_subject.CHROMS:
                continue
            for vcf_record in vcf_subject.CHROMS[chroms_index]:

                interval = self.assign_to_window(intervals, vcf_record)
                if not interval:
                    continue
                intervals[interval].append(vcf_record)

            self.pms_dict[chroms_index] = {}
            self.pms_burden[chroms_index] = {}
            
            for interval in intervals:
                ms_tmp, bad_records_tmp = self.compute_ms(ref_chrom.seq, intervals[interval])
                self.bad_records += bad_records_tmp
                if len(intervals[interval]) >= 96:
                    self.pms_dict[chroms_index][(interval.start, interval.end)] = ms_tmp
                self.pms_burden[chroms_index][(interval.start, interval.end)] = len(intervals[interval])
                
                
                
    def window_generator(self, chrom_length, window_size):
    
        intervals = {}
        for i in range(1, chrom_length + 1, window_size):
            intervals[Interval(i, i + window_size - 1)] = []
        if i + window_size - 1 < chrom_length + 1:
            intervals[Interval(i + window_size, chrom_length)] = []
        return intervals
    
    def assign_to_window(self, intervals, vcf_record):
        for interval in intervals:
            if interval.is_included(vcf_record[0]):
                return interval
        return None
    
    def compute_ms(self, SEQ, Mutations):
        
        MATE_PAIR = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        MS = {'ACA-A':0,'ACC-A':0,'ACG-A':0,'ACT-A':0,'CCA-A':0,
              'CCC-A':0,'CCG-A':0,'CCT-A':0,'GCA-A':0,'GCC-A':0,
              'GCG-A':0,'GCT-A':0,'TCA-A':0,'TCC-A':0,'TCG-A':0,
              'TCT-A':0,'ACA-G':0,'ACC-G':0,'ACG-G':0,'ACT-G':0,
              'CCA-G':0,'CCC-G':0,'CCG-G':0,'CCT-G':0,'GCA-G':0,
              'GCC-G':0,'GCG-G':0,'GCT-G':0,'TCA-G':0,'TCC-G':0,
              'TCG-G':0,'TCT-G':0,'ACA-T':0,'ACC-T':0,'ACG-T':0,
              'ACT-T':0,'CCA-T':0,'CCC-T':0,'CCG-T':0,'CCT-T':0,
              'GCA-T':0,'GCC-T':0,'GCG-T':0,'GCT-T':0,'TCA-T':0,
              'TCC-T':0,'TCG-T':0,'TCT-T':0,'ATA-A':0,'ATC-A':0,
              'ATG-A':0,'ATT-A':0,'CTA-A':0,'CTC-A':0,'CTG-A':0,
              'CTT-A':0,'GTA-A':0,'GTC-A':0,'GTG-A':0,'GTT-A':0,
              'TTA-A':0,'TTC-A':0,'TTG-A':0,'TTT-A':0,'ATA-C':0,
              'ATC-C':0,'ATG-C':0,'ATT-C':0,'CTA-C':0,'CTC-C':0,
              'CTG-C':0,'CTT-C':0,'GTA-C':0,'GTC-C':0,'GTG-C':0,
              'GTT-C':0,'TTA-C':0,'TTC-C':0,'TTG-C':0,'TTT-C':0,
              'ATA-G':0,'ATC-G':0,'ATG-G':0,'ATT-G':0,'CTA-G':0,
              'CTC-G':0,'CTG-G':0,'CTT-G':0,'GTA-G':0,'GTC-G':0,
              'GTG-G':0,'GTT-G':0,'TTA-G':0,'TTC-G':0,'TTG-G':0,'TTT-G':0}
        # Include all of the mutations
        bad_records_tmp = 0
        for ii in range(len(Mutations)):
            vcf_record = Mutations[ii]
            POS = vcf_record[0] - 1
            REF = vcf_record[1]
            ALT = vcf_record[2]

            if POS >= len(SEQ) or POS < 0:
                bad_records_tmp += 1
                continue
            # Error check
            if SEQ[POS] != REF:
                bad_records_tmp += 1
                continue
    #             raise ValueError('Error: Reference inconsistent')

            # Generate dictionary key
            KEY = SEQ[(POS-1):(POS+2)] + '-'  + ALT

            # Change key for A or G
            if (SEQ[POS]=='A') or (SEQ[POS]=='G'):
                if self.mode == 0:
                    bad_records_tmp += 1
                    continue
                KEY = MATE_PAIR[KEY[0]] + MATE_PAIR[KEY[1]] + MATE_PAIR[KEY[2]] + '-' + MATE_PAIR[KEY[4]]
            MS[KEY] += 1
            self.overall_mutational_signature[KEY] += 1
        
        return MS, bad_records_tmp
    
    def organize_ms_burden_mat(self, ms_burden, WINDOW_SIZE = 10**6):
    
        chrom_len = max([len(ms_burden[chrom]) for chrom in ms_burden])
        ms_burden_mat = -1 * np.ones((22, chrom_len))

        for chrom_idx in ms_burden:
            for start, end in ms_burden[chrom_idx]:
                col_idx = end // WINDOW_SIZE
                ms_burden_mat[chrom_idx - 1, col_idx - 1] = ms_burden[chrom_idx][(start, end)]
        return ms_burden_mat

    def plot_ms_burden(self, ms_burden_mat, subject_idx = 1):

        masked_array = np.ma.array(ms_burden_mat, mask=(ms_burden_mat == -1))
        fig = plt.figure(figsize = (40, 100))
        cmap = matplotlib.cm.hot
        cmap.set_bad('white',1.)
        plt.imshow(masked_array, interpolation='nearest', cmap=cmap)
        plt.xticks(np.arange(0, ms_burden_mat.shape[1], 50), [str(i) + 'mbps' for i in range(0, 250, 50)])
        plt.yticks(np.arange(22), ['ch' + str(i) for i in range(1, 23)])
        plt.title('Mutational Burden subj_' + str(subject_idx))
        fig.savefig('../Figure/mutational_burden_subj_' + str(subject_idx) + '.pdf', format='pdf', dpi=300, bbox_inches='tight')
        plt.close('all')
        # plt.show()
        
    def visualize_ms_burden(self, subject_id):
        ms_burden_mat = self.organize_ms_burden_mat(self.pms_burden)
        self.plot_ms_burden(ms_burden_mat, subject_id)
        
        
    def visualize_ms_plot(self, MS, sub_title):
        fig = plt.figure(figsize = (15, 8))
        ms_keys = sorted(list(MS.keys()), key= lambda x: (x[1], x[4], x[0], x[2]))
        ms_values = [MS[key] for key in ms_keys]
        barlist = plt.bar(range(len(MS.values())),ms_values,color=[0.5,0.5,0.8])
        for ii in range(0,16):
            barlist[ii].set_color(np.array([135,206,250])/256)
        for ii in range(16,16*2):
            barlist[ii].set_color(np.array([0,0,0]))
        for ii in range(16*2,16*3):
            barlist[ii].set_color(np.array([1, 0, 0]))
        for ii in range(16*3,16*4):
            barlist[ii].set_color(np.array([0.5, 0.5, 0.5]))
        for ii in range(16*4,16*5):
            barlist[ii].set_color(np.array([0,1,0]))
        for ii in range(16*5,16*6):
            barlist[ii].set_color(np.array([255,182,193])/256)
        plt.ylabel('Frequency')
        plt.title('Mutation Signature ' + sub_title)
        tick_labels = ms_keys
        plt.xticks(range(96),tick_labels,rotation=90,fontsize=7)
        fig.savefig('../Figure/mutational_signature_' + sub_title + '.pdf', format='pdf', dpi=300, bbox_inches='tight')
        plt.close('all')
        # plt.show()
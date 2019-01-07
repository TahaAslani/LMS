def main():
    from Positional_Mutational_Signature import Positional_Mutational_Signature
    for subject_id in range(1, 8):
        print('Computing PMS for subject %d...' % subject_id)
        pms = Positional_Mutational_Signature(subject_file_path = '../Texas_vcfs/wheeljack' + str(subject_id) + '.vcf', mode = '0')
        pms.visualize_ms_burden(subject_id)
        pms.visualize_ms_plot(pms.overall_mutational_signature, 'subj_' + str(subject_id) + '_overall')
        for chrom_idx in pms.pms_dict:
            for interval in pms.pms_dict[chrom_idx]:
                sub_title = 'subj_%d_chr_%d_from_%d_to_%d' % (subject_id, chrom_idx, interval[0], interval[1])
                pms.visualize_ms_plot(pms.pms_dict[chrom_idx][interval], sub_title)
        print('bad_reads', pms.bad_records)

if __name__ == '__main__':
    main()
''' '''

from unittest import TestCase
import gwass


class TestSummaryStatistics(TestCase):

    def setUp(self):
        ''' '''
        summary_statistics_path = 'gwass/tests/test_data/summary_statistics/summary_statistics_test.tsv' 
        snps_path = 'gwass/tests/test_data/summary_statistics/snps_test.tsv'
        self.ss = gwass.SummaryStatistics(snps_path=snps_path, summary_statistics_path=summary_statistics_path,
                                          snp_col='SNP', a1_col='A1', a2_col='A2', 
                                          beta_hat_col='BETA', se_col='SE', p_value_col='NA', 
                                          sample_size_col='NA', regression_type_col='linear')

    def test_read_summary_statistics(self):
        ''' '''
        self.ss._read_summary_statistics()
        self.assertEqual(self.ss.summary_statistics.shape[0], 6) 

    def test_read_snps(self):
        ''' '''
        self.ss._read_snps()
        self.assertEqual(self.ss.snps.shape[0], 4)

    def test_alleles_to_upper(self):
        ''' '''
        self.ss._read_summary_statistics()
        self.ss._alleles_to_upper()
        self.assertEqual(self.ss.summary_statistics['a1'].tolist()[0], 'C')        
        self.assertEqual(self.ss.summary_statistics['a2'].tolist()[0], 'G')    

    def test_remove_indels(self):
        ''' '''
        self.ss._read_summary_statistics()
        self.ss._alleles_to_upper()
        self.ss._remove_indels()
        self.assertEqual(self.ss.summary_statistics.shape[0], 5) 

    def test_merge_summary_statistics_with_snps(self):
        ''' '''
        self.ss._read_summary_statistics()
        self.ss._read_snps()
        self.ss._alleles_to_upper()
        self.ss._remove_indels()
        self.ss._merge_summary_statistics_with_snps()
        self.assertEqual(self.ss.summary_statistics.shape[0], 4) 
        self.assertEqual(self.ss.summary_statistics.shape[1], 25) 

    def test_remove_strand_ambiguous_snp(self):
        ''' '''
        self.ss._read_summary_statistics()
        self.ss._read_snps()
        self.ss._alleles_to_upper()
        self.ss._remove_indels()
        self.ss._merge_summary_statistics_with_snps()
        self.ss._remove_strand_ambiguous_snps()
        self.assertEqual(self.ss.summary_statistics.shape[0], 2) 
         
    def test_orient_snp_effect_sign(self):
        ''' '''
        self.ss._read_summary_statistics()
        self.ss._read_snps()
        self.ss._alleles_to_upper()
        self.ss._remove_indels()
        self.ss._merge_summary_statistics_with_snps()
        self.ss._remove_strand_ambiguous_snps()
        self.ss._orient_snp_effect_signs()
        self.assertEqual(self.ss.summary_statistics['beta_hat'].tolist()[0], .1)

    def test_clean_data_frame(self):
        ''' '''
        self.ss._read_summary_statistics()
        self.ss._read_snps()
        self.ss._alleles_to_upper()
        self.ss._remove_indels()
        self.ss._merge_summary_statistics_with_snps()
        self.ss._remove_strand_ambiguous_snps()
        self.ss._orient_snp_effect_signs()
        self.ss._clean_data_frame()
        header = list(self.ss.summary_statistics.columns.values)
        self.assertEqual('a1' in header, False)


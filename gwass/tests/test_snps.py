from unittest import TestCase

import gwass
import pysam


class TestSnp(TestCase):

    def setUp(self):
        ''' ''' 
        # variant files
        snp_aa_vcf = pysam.VariantFile('gwass/tests/test_data/snp_aa.vcf.gz')    
        snp_no_aa_0_vcf = pysam.VariantFile('gwass/tests/test_data/snp_no_aa_0.vcf.gz')
        snp_no_aa_1_vcf = pysam.VariantFile('gwass/tests/test_data/snp_no_aa_1.vcf.gz')
        snp_small_aa_vcf = pysam.VariantFile('gwass/tests/test_data/snp_small_aa.vcf.gz')
        # records
        snp_aa_rec = [rec for rec in snp_aa_vcf][0]
        snp_no_aa_0_rec = [rec for rec in snp_no_aa_0_vcf][0]
        snp_no_aa_1_rec = [rec for rec in snp_no_aa_1_vcf][0]
        snp_small_aa_rec = [rec for rec in snp_small_aa_vcf][0]
        # snps
        self.snp_aa = gwass.Snp(snp_aa_rec)
        self.snp_no_aa_0 = gwass.Snp(snp_no_aa_0_rec)
        self.snp_no_aa_1 = gwass.Snp(snp_no_aa_1_rec)
        self.snp_small_aa = gwass.Snp(snp_small_aa_rec)

    def test_coordinates(self):
        ''' ''' 
        # snp_aa tests 
        self.assertEqual(self.snp_aa.chrom, '1')
        self.assertEqual(self.snp_aa.pos, 524937) 
        self.assertEqual(self.snp_aa.rsid, 'rs556300838')
        # snp_no_aa_0 tests 
        self.assertEqual(self.snp_no_aa_0.chrom, '1')
        self.assertEqual(self.snp_no_aa_0.pos, 91264) 
        self.assertEqual(self.snp_no_aa_0.rsid, 'rs550687227')
        # snp_no_aa_1         
        self.assertEqual(self.snp_no_aa_1.chrom, '1')
        self.assertEqual(self.snp_no_aa_1.pos, 54727) 
        self.assertEqual(self.snp_no_aa_1.rsid, 'rs531848864')
        # snp_small_aa 
        self.assertEqual(self.snp_small_aa.chrom, '1')
        self.assertEqual(self.snp_small_aa.pos, 52687) 
        self.assertEqual(self.snp_small_aa.rsid, 'rs550953937')
    
    def test_ref_alt_alleles(self):
        ''' '''
        # snp_aa
        self.assertEqual(self.snp_aa.ref_allele, 'C')
        self.assertEqual(self.snp_aa.alt_allele, 'T')
        # snp_no_aa_0
        self.assertEqual(self.snp_no_aa_0.ref_allele, 'A')
        self.assertEqual(self.snp_no_aa_0.alt_allele, 'C')
        # snp_no_aa_1
        self.assertEqual(self.snp_no_aa_1.ref_allele, 'T')
        self.assertEqual(self.snp_no_aa_1.alt_allele, 'C')
        # snp_small_aa
        self.assertEqual(self.snp_small_aa.ref_allele, 'T')
        self.assertEqual(self.snp_small_aa.alt_allele, 'C')

    def test_allele_type(self):
        ''' '''
        # snp_aa
        self.assertEqual(self.snp_aa.allele_type, 'derived_ancestral')
        # snp_no_aa_0
        self.assertEqual(self.snp_no_aa_0.allele_type, 'alternate_reference')
        # snp_no_aa_1
        self.assertEqual(self.snp_no_aa_1.allele_type, 'alternate_reference')
        # snp_small_aa
        self.assertEqual(self.snp_small_aa.allele_type, 'derived_ancestral')

    def test_get_ancestral_derived_alleles(self):
        ''' '''
        # snp_aa 
        self.assertEqual(self.snp_aa.ancestral_allele, 'C')
        self.assertEqual(self.snp_aa.derived_allele, 'T')    
        # snp_no_aa_0  
        self.assertEqual(self.snp_no_aa_0.ancestral_allele, 'NA')
        self.assertEqual(self.snp_no_aa_0.derived_allele, 'NA')    
        # snp_no_aa_1 
        self.assertEqual(self.snp_no_aa_1.ancestral_allele, 'NA')
        self.assertEqual(self.snp_no_aa_1.derived_allele, 'NA')    
        # snp_small_aa 
        self.assertEqual(self.snp_small_aa.ancestral_allele, 'T')
        self.assertEqual(self.snp_small_aa.derived_allele, 'C')    

    def test_get_effect_other_alleles(self):
        ''' '''
        # snp_aa
        self.assertEqual(self.snp_aa.effect_allele, 'T')
        self.assertEqual(self.snp_aa.other_allele, 'C')
        # snp_no_aa_0
        self.assertEqual(self.snp_no_aa_0.effect_allele, 'C')
        self.assertEqual(self.snp_no_aa_0.other_allele, 'A')
        # snp_no_aa_1
        self.assertEqual(self.snp_no_aa_1.effect_allele, 'C')
        self.assertEqual(self.snp_no_aa_1.other_allele, 'T')
        # snp_small_aa
        self.assertEqual(self.snp_small_aa.effect_allele, 'C')
        self.assertEqual(self.snp_small_aa.other_allele, 'T')

    def test_get_allele_frequencies(self):
        ''' '''
        # snp_aa
        self.assertAlmostEqual(self.snp_aa.f_sas, 0.0072, places=4)  
        self.assertAlmostEqual(self.snp_aa.f_afr, 0.0, places=4)  
        self.assertAlmostEqual(self.snp_aa.f_eur, 0.0, places=4)  
        self.assertAlmostEqual(self.snp_aa.f_eas, 0.0, places=4)  
        self.assertAlmostEqual(self.snp_aa.f_amr, 0.0, places=4)  
        # snp_no_aa_0
        self.assertAlmostEqual(self.snp_no_aa_0.f_sas, 0.0072, places=4)  
        self.assertAlmostEqual(self.snp_no_aa_0.f_afr, 0.0, places=4)  
        self.assertAlmostEqual(self.snp_no_aa_0.f_eur, 0.0089, places=4)  
        self.assertAlmostEqual(self.snp_no_aa_0.f_eas, 0.0, places=4)  
        self.assertAlmostEqual(self.snp_no_aa_0.f_amr, 0.0115, places=4)  
        # snp_no_aa_1
        self.assertAlmostEqual(self.snp_no_aa_1.f_sas, 0.001, places=4)  
        self.assertAlmostEqual(self.snp_no_aa_1.f_afr, 0.0, places=4)  
        self.assertAlmostEqual(self.snp_no_aa_1.f_eur, 0.0, places=4)  
        self.assertAlmostEqual(self.snp_no_aa_1.f_eas, 0.0, places=4)  
        self.assertAlmostEqual(self.snp_no_aa_1.f_amr, 0.0, places=4)  
        # snp_small_aa        
        self.assertAlmostEqual(self.snp_small_aa.f_sas, 0.0, places=4)  
        self.assertAlmostEqual(self.snp_small_aa.f_afr, 0.0008, places=4)  
        self.assertAlmostEqual(self.snp_small_aa.f_eur, 0.0, places=4)  
        self.assertAlmostEqual(self.snp_small_aa.f_eas, 0.0, places=4)  
        self.assertAlmostEqual(self.snp_small_aa.f_amr, 0.0, places=4)  



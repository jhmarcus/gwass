''' TODO: description '''

import sys
import pysam


class Snp:
    '''
    extract ancestral and derived allele information from
    1kg_phase3 sites file for a snp based on pysam record

    TODO: add neighboring ancestral alleles for mutation context analyses
    '''

    def __init__(self, rec, ref_genome):

        # pysam record
        self.rec = rec

        # pysam fasta
        self.ref_genome = ref_genome

        # snp attributes
        self.chrom = self.rec.chrom
        self.pos = self.rec.pos
        self.rsid = self.rec.id
        self.ref_allele = self.rec.ref
        self.alt_allele = self.rec.alts[0]
        self.minor_allele = None
        self.major_allele = None
        self.ancestral_allele = None
        self.derived_allele = None

        # allele frequencies
        self.f_sas = None
        self.f_afr = None
        self.f_eas = None
        self.f_amr = None
        self.f_eur = None

        # mut context
        self.ref_l2 = None
        self.ref_l1 = None
        self.ref_r1 = None
        self.ref_r2 = None

        # bases
        self.bases = ['A', 'C', 'T', 'G']

        # set attributes
        self._get_minor_major_alleles()
        self._get_ancestral_derived_alleles()
        self._get_allele_frequencies()
        self._get_mutation_context()

    def _get_minor_major_alleles(self):
        ''' TODO: test '''
        if self.rec.info['AF'][0] < .5:
            self.minor_allele = self.alt_allele
            self.major_allele = self.ref_allele
        else:
            self.minor_allele = self.ref_allele
            self.major_allele = self.alt_allele

    def _get_ancestral_derived_alleles(self):
        ''' extract ancestral and derived allele from vcf annotation '''

        # check if ancestral allele is in record
        if 'AA' in list(self.rec.info):
            self.ancestral_allele = self.rec.info['AA'][0].upper()

            # check if ancestral allelle is a valid base / is the ref or alt base double mutation?
            if ((self.ancestral_allele in self.bases) and
                    (self.ancestral_allele in [self.ref_allele, self.alt_allele])):

                # check if ancestral allele is reference allele
                if self.ref_allele == self.ancestral_allele:
                    self.derived_allele = self.alt_allele
                else:
                    self.derived_allele = self.ref_allele
            else:
                self.ancestral_allele = 'NA'
                self.derived_allele = 'NA'
        else:
            self.ancestral_allele = 'NA'
            self.derived_allele = 'NA'

    def _get_allele_frequencies(self):
        '''
        '''
        self.f_sas = 1.0 - self.rec.info['SAS_AF'][0]
        self.f_afr = 1.0 - self.rec.info['AFR_AF'][0]
        self.f_eas = 1.0 - self.rec.info['EAS_AF'][0]
        self.f_eur = 1.0 - self.rec.info['EUR_AF'][0]
        self.f_amr = 1.0 - self.rec.info['AMR_AF'][0]

    def _get_mutation_context(self):
        '''
        get references bases flanking the snp for mutation context

        TODO: test!
        '''
        mut_context = self.ref_genome.fetch(self.chrom, self.pos-3, self.pos+2)
        self.ref_l2 = mut_context[0]
        self.ref_l1 = mut_context[1]
        self.ref_r1 = mut_context[3]
        self.ref_r2 = mut_context[4]


class Snps:
    '''
    extract ancestral and derived allele information from 1kg_phase3
    sites file for all snps into a pandas DataFrame
    '''

    def __init__(self, sites_vcf_path, ref_genome_path):

        # input
        self.sites_vcf_path = sites_vcf_path
        self.ref_genome_path = ref_genome_path
        self.snps = None

    def write_snps(self):
        '''writes snps to file'''

        vcf = pysam.VariantFile(self.sites_vcf_path)
        ref_genome = pysam.FastaFile(self.ref_genome_path)
        header = 'chrom\tpos\tsnp\tref_allele\talt_allele\tmajor_allele\tminor_allele\t' \
                 'ancestral_allele\tderived_allele\tref_base_l2\tref_base_l1\tref_base_r1\tref_base_r2\t' \
                 'f_sas\tf_afr\tf_eas\tf_eur\tf_amr\n'
        sys.stdout.write(header)
        for rec in vcf.fetch():
            snp = Snp(rec, ref_genome)
            sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t' \
                             '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(snp.chrom, snp.pos, snp.rsid,
                                                                           snp.ref_allele, snp.alt_allele,
                                                                           snp.major_allele, snp.minor_allele,
                                                                           snp.ancestral_allele, snp.derived_allele,
                                                                           snp.ref_l2, snp.ref_l1, snp.ref_r1, snp.ref_r2,
                                                                           float('%.5f' % snp.f_sas), float('%.5f' % snp.f_afr),
                                                                           float('%.5f' % snp.f_eas), float('%.5f' % snp.f_eur),
                                                                           float('%.5f' % snp.f_amr)))

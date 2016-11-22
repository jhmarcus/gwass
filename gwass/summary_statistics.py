''' TODO: description '''

import pandas as pd
import numpy as np
import itertools as it

class SummaryStatistics: 
    '''class for processing genome-wide association summary statistics'''

    def __init__(self, summary_statistics_path, snps_path):
        
        # input 
        self.summary_statistics_path = summary_statistics_path
        self.snps_path = snps_path

        # dict helpers
        self.base_complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        self.bases = self.base_complement.keys()
        self.strand_ambiguous_alleles = {''.join(x): x[0] == self.base_complement[x[1]]
                                         for x in it.product(self.bases, self.bases)
                                         if x[0] != x[1]} 
        self.default_col_names = {
                                  'snp': 'snp', # snp
                                  'SNP':  'snp', 
                                  'snpid': 'snp', 
                                  'SNPID': 'snp',
                                  'MarkerName': 'snp',
                                  'rs_number': 'snp',
                                  'a1': 'a1', # a1
                                  'reference': 'a1',
                                  'reference_allele': 'a1',
                                  'A1': 'a1', 
                                  'Allele1': 'a1', 
                                  'allele1': 'a1',
                                  'effect_allele': 'a1',
                                  'a2': 'a2', # a2
                                  'other_allele': 'a2',
                                  'A2': 'a2',
                                  'Allele2': 'a2', 
                                  'allele2': 'a2',
                                  'beta': 'beta_hat', # beta
                                  'BETA': 'beta_hat', 
                                  'b': 'beta_hat',
                                  'Effect': 'beta_hat',
                                  'effect': 'beta_hat',
                                  'or': 'or', # or
                                  'OR': 'or',
                                  'se': 'se', # se
                                  'SE': 'se',
                                  'StdErr': 'se',
                                  'stderr': 'se' 
                                } 


        # DataFrames
        self.summary_statistics = None
        self.snps = None
        
        # fil counts
        self.n_snps = None
        self.n_fil_indels = None
        self.n_fil_merged = None
        self.n_fil_strand = None
        self.n_fil_snps = None

    def _read_summary_statistics(self):
        '''reads summary statistics and converts column names to default schema'''
        # read header and table
        df_header = list(pd.read_table(self.summary_statistics_path, nrows=0, index_col=False, sep=None, engine='python'))
        df = pd.read_table(self.summary_statistics_path, sep='\t', index_col=False, skiprows=1, header=None)  
        df.columns = df_header 
        df_header = [col_name for col_name in df_header if col_name in self.default_col_names.keys()] 
        if len(df_header) == 0:
            raise ValueError('summary statistics file is empty or none of column names match default column names.')
        df = df[df_header] 

        # convert header to default schema
        df_header = [self.default_col_names[col_name] for col_name in list(df)] 
        df.columns = df_header
        if ('or' not in df.columns) and ('beta_hat' not in df.columns):
            raise ValueError('Missing or / beta_hat column from summary statistics file or or / beta_hat column does not match default column names')
        for col in ['snp', 'a1', 'a2', 'se']:
            if col not in df.columns:
                raise ValueError('Missing {} from summary statistics file or snp column does not match default column names'.format(col)) 

        # assign to df summary_statistics
        self.summary_statistics = df 

        # snp counts 
        self.n_snps = self.summary_statistics.shape[0] 
        self.n_fil_snps = self.summary_statistics.shape[0]

    def _read_snps(self):
        '''reads snps'''
        self.snps = pd.read_table(self.snps, sep='\t', index_col=False)
            
    def _alleles_to_upper(self):
        '''converts alleles to upper case'''
        self.summary_statistics['a1'] = self.summary_statistics['a1'].apply(lambda x: x.upper())
        self.summary_statistics['a2'] = self.summary_statistics['a2'].apply(lambda x: x.upper())

    def _remove_indels(self):
        '''remove indels from summary statistics'''
        self.summary_statistics = self.summary_statistics.loc[self.summary_statistics['a1'].isin(self.bases)]
        self.summary_statistics = self.summary_statistics.loc[self.summary_statistics['a2'].isin(self.bases)] 
        
        # snp counts
        self.n_fil_indels = self.n_fil_snps - self.summary_statistics.shape[0]
        self.n_fil_snps = self.summary_statistics.shape[0]

    def _merge_summary_statistics_with_snps(self):
        '''joins summary statistics and snps'''
        self.summary_statistics = self.summary_statistics.merge(self.snps, on='snp', how='inner') 
        
        # snp counts
        self.n_fil_merged = self.n_fil_snps - self.summary_statistics.shape[0]
        self.n_fil_snps = self.summary_statistics.shape[0]        

    def _transform_or_to_beta_hat(self):
        '''log transform of odds ratios to convert to beta_hat'''
        self.summary_statistics['beta'] = np.log(self.summary_statistics['or']) 
        self.summary_statistics.drop('or', axis=1, inplace=True)
        
    def _remove_strand_ambiguous_snp(a1, a2, strand_ambiguous_alleles):
        ''' '''
        a1a2 = a1 + a2
        if strand_ambiguous_alleles[a1a2]:
            return(False)
        else:
            return(True)
    
    def _remove_strand_ambiguous_snps(self):
        ''' '''
        self.summary_statistics = self.summary_statistics[self.summary_statistics.apply(lambda row: self.remove_strand_ambiguous_snp(row['a1'], 
                                                                                        row['a2'], self.strand_ambiguous_alleles),  axis=1)]
        # snp counts
        self.n_fil_strand = self.n_fil_snps - self.summary_statistics.shape[0]        
        self.n_fil_strand = self.summary_statistics.shape[0]

    '''
    TODO: add rest of code and function run all qc steps
    '''
            


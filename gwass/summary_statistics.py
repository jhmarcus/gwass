''' TODO: description '''

import pandas as pd
import numpy as np
import itertools as it


class SummaryStatistics: 
    '''class for processing genome-wide association summary statistics'''

    def __init__(self, summary_statistics_path, snps_path, snp_col, 
                 a1_col, a2_col, beta_hat_col, se_col, p_value_col, 
                 sample_size_col, regression_type_col):
        
        # input 
        self.summary_statistics_path = summary_statistics_path
        self.snps_path = snps_path

        # dict helpers
        self.base_complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        self.bases = self.base_complement.keys()
        self.strand_ambiguous_alleles = {''.join(x): x[0] == self.base_complement[x[1]]
                                         for x in it.product(self.bases, self.bases)
                                         if x[0] != x[1]} 

        # default colnames
        self.regression_type_col = regression_type_col
        self.p_value_col = p_value_col
        self.sample_size_col = sample_size_col
        self.default_col_names = {
                                  snp_col: 'snp',
                                  a1_col: 'a1',
                                  a2_col: 'a2',
                                  beta_hat_col: 'beta_hat',
                                  se_col: 'se',
                                  p_value_col: 'p_value',
                                  sample_size_col: 'sample_size'
                                 }
 
        # DataFrames
        self.summary_statistics = None
        self.snps = None
        
        # fil counts
        self.n_snps = None
        self.n_fil_indels = None
        self.n_fil_merged = None
        self.n_fil_strand = None
        self.n_fil_inconsistent = None
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

        # colname assertions
        if ('or' not in df.columns) and ('beta_hat' not in df.columns):
            raise ValueError('Missing or / beta_hat column from summary statistics file or or / beta_hat column does not match default column names')
        for col in ['snp', 'a1', 'a2', 'se']:
            if col not in df.columns:
                raise ValueError('Missing {} from summary statistics file or snp column does not match default column names'.format(col)) 
        
        # add p_value and sample size col if missing 
        if self.p_value_col == 'NA':
          df['p_value'] = 'NA'
        if self.sample_size_col == 'NA':
          df['sample_size'] = 'NA'

        # assign to df summary_statistics
        self.summary_statistics = df 

        # snp counts 
        self.n_snps = self.summary_statistics.shape[0] 
        self.n_fil_snps = self.summary_statistics.shape[0]

    def _read_snps(self):
        '''reads snps'''

        self.snps = pd.read_table(self.snps_path, sep='\t', index_col=False)
            
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

        self.summary_statistics['beta_hat'] = np.log(self.summary_statistics['beta_hat']) 
    
    def _remove_strand_ambiguous_snps(self):
        '''removes snps where strand cannot be determined'''

        self.summary_statistics = self.summary_statistics[self.summary_statistics.apply(lambda row: remove_strand_ambiguous_snp(row['a1'], row['a2'], 
                                                                                                                                self.strand_ambiguous_alleles), 
                                                                                                                                axis=1)]
        # snp counts
        self.n_fil_strand = self.n_fil_snps - self.summary_statistics.shape[0]        
        self.n_fil_snps = self.summary_statistics.shape[0]

    def _remove_inconsistent_snps(self):
        '''removes snps where alleles are inconsistent between both datasets'''

        self.summary_statistics = self.summary_statistics[self.summary_statistics.apply(lambda row: remove_inconsistent_snp(row['a1'], row['a2'], 
                                                                                                                            row['ref_allele'], row['alt_allele'],
                                                                                                                            self.base_complement), axis=1)]
        # snp counts
        self.n_fil_inconsistent = self.n_fil_snps - self.summary_statistics.shape[0] 
        self.n_fil_snps = self.summary_statistics.shape[0]

    def _orient_snp_effect_signs(self):
        '''applys _orient_snp_effect_sign to all snps in the DataFrame'''

        self.summary_statistics['beta_hat'] = self.summary_statistics.apply(lambda row: orient_snp_effect_sign(row['snp'], row['a1'], row['a2'], 
                                                                                                               row['ref_allele'], row['alt_allele'], 
                                                                                                               row['beta_hat'], self.base_complement), axis=1)

    def _clean_data_frame(self):
        '''drops uncessary cols, reorders cols and sorts by chrom-pos'''

        # drop a1 and a2 columns
        self.summary_statistics.drop('a1', axis=1, inplace=True) 
        self.summary_statistics.drop('a2', axis=1, inplace=True) 

        # sort columns
        self.summary_statistics = self.summary_statistics.sort_values(['chrom', 'pos'], ascending=[True, True])

    def clean_summary_statistics(self, out):
        '''public call of functions to clean summary statistics'''

        # summary_statistics pipleline
        self._read_summary_statistics()
        self._read_snps()
        self._alleles_to_upper()
        self._remove_indels()
        self._merge_summary_statistics_with_snps()
        if self.regression_type_col == 'logistic':
            self._transform_or_to_beta_hat()
        self._remove_strand_ambiguous_snps()
        self._remove_inconsistent_snps()
        self._orient_snp_effect_signs()
        self._clean_data_frame()
        col_names = ['chrom', 'pos', 'snp', 'ref_allele', 'alt_allele', 'beta_hat', 'se', 'p_value', 'sample_size', 'f_eur']
        self.summary_statistics[col_names].to_csv('{}.tsv.gz'.format(out), sep='\t', compression='gzip', index=False, na_rep='NA', columns=col_names)

        # snp counts 
        fil_df = pd.DataFrame({'n_snps': np.array([self.n_snps]), 'n_fil_snps': np.array([self.n_fil_snps]), 'n_fil_indels': np.array([self.n_fil_indels]), 
                               'n_fil_merged': np.array([self.n_fil_merged]), 'n_fil_strand': np.array([self.n_fil_strand]), 
                               'n_fil_inconsistent': np.array([self.n_fil_inconsistent])})

        fil_df_cols = ['n_snps', 'n_fil_snps', 'n_fil_indels', 'n_fil_merged', 'n_fil_strand', 'n_fil_inconsistent']
        fil_df.to_csv('{}_fil.tsv'.format(out), sep='\t', index=False, na_rep='NA', columns=fil_df_cols)

# helper functions for methods

def remove_strand_ambiguous_snp(a1, a2, strand_ambiguous_alleles):
    '''remove snps where strand cannot be determined'''

    a1a2 = a1 + a2
    if strand_ambiguous_alleles[a1a2]:
        return(False)
    else:
        return(True)

def remove_inconsistent_snp(a1, a2, effect_allele, other_allele, base_complement):
    '''remove snps with inconsistent alleles between the datasets'''

    if (a1 == other_allele) and (a2 == effect_allele):    
        return(True)
    elif (a1 == base_complement[other_allele]) and (a2 == base_complement[effect_allele]):
        return(True)
    elif (a1 == effect_allele) and (a2 == other_allele):
        return(True)
    elif (a1 == base_complement[effect_allele]) and (a2 == base_complement[other_allele]):
        return(True)
    else:
        return(False)

def orient_snp_effect_sign(snp, a1, a2, effect_allele, other_allele, beta_hat, base_complement):
    '''
    orients beta_hat to measure effect on derived allele 
    when known otherwise to the alternate allele.
    '''

    if (a1 == other_allele) and (a2 == effect_allele):    
        beta_hat = -beta_hat
    elif (a1 == base_complement[other_allele]) and (a2 == base_complement[effect_allele]):
        beta_hat = -beta_hat
    elif (a1 == effect_allele) and (a2 == other_allele):
        beta_hat = beta_hat
    elif (a1 == base_complement[effect_allele]) and (a2 == base_complement[other_allele]):
        beta_hat = beta_hat
    else:
        raise ValueError('inconsistent alleles: snp={},a1={},a2={},effect_allele={},other_allele={}'.format(snp, a1, a2, effect_allele, other_allele))

    return(beta_hat)


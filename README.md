# gwass

*in progress*  

A python package for analyzing Genome-Wide Association Summary Statistics (gwass). Currently, the package has two primary functions:

* Creating a table of SNPs with meta-data about derived and ancestral alleles / frequencies  
* Cleaning up GWAS summary statistics and orienting the sign of effects to derived alleles

Currently the package has a specific focus such that it only works with the 1000 Genomes Project data but in the future this will be generalized.

## Pre-requitesites

* `python3`
* `pysam`
* `numpy`
* `pandas`

## Installation

I recommend creating and activating a new conda python environment:

```
conda create -n <env_name> python
source activate <env_name>
```

Clone the repository and install into your python env

```
git clone https://github.com/jhmarcus/gwass
cd gwass/
make install
```

Make sure unit tests are working:

```
make test
```

## Demo

### Download 1KG sites vcf and reference genome

```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz
```

### Create sites vcf with only bi-allelic SNPs

Assumes bcftools is installed (TODO: add bcftools command)

```
bcftools
```

#### Download Height GWAS summary statistics

```
wget http://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz
```

### Create SNPs

Once the package has been installed it should add command-line executables to your path. One of these scripts will create a SNPs table given the 1KG sites vcf:

```
create-snps --sites <path>/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.biallelic.vcf.gz \
            --ref <path>/human_g1k_v37.fasta | \
            gzip -c > 1kg_phase3_snps.tsv.gz
```

*note that this will take some time and is a pretty large file*. The script will create the file `1kg_phase3_snps.tsv.gz` with columns:

* `chrom`: chromsome number
* `pos`: reference genome position
* `snp`: rsID
* `effect_allele`: allele which effect sign is to be oriented to
* `other_allele`: other allele which effect sign is not oriented to
* `derived_allele`: the derived allele determined by multi-way alignment
* `ancestral_allele` the ancestral allele determined by multi-way alignment
* `minor_allele`: the globally minor allele determined by 1KG Phase3
* `major_allele`: the globally major allele determined by 1KG Phase3
* `ref_allele`: the allele in the reference genome
* `alt_allele`: the alternate allele to the reference genome
* `allele_type`: the type of the effect allele / other allele which is either derived_ancestral or minor_major
* `ref_base_l2`: the reference base two positions to the left of the SNP
* `ref_base_l1`: the reference base one positions to the left of the SNP
* `ref_base_r1`: the reference base one positions to the right of the SNP
* `ref_base_r2`: the reference base two positions to the right of the SNP
* `f_sas`: effect allele frequency in 1KG Phase3 south asians
* `f_eas`: effect allele frequency in 1KG Phase3 east asians
* `f_eur`: effect allele frequency in 1KG Phase3 europeans
* `f_afr`: effect allele frequency in 1KG Phase3 africans
* `f_amr`: effect allele frequency in 1KG Phase3 americans

### Clean summary statistics

*note that this will take a bit of memory if loading the whole SNPs file*

```
clean-summary-statistics --summary_statistics <path>/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz \
                         --snps <path>/1kg_phase3_snps  
                         --out giant_height_summary_statistics
```

This will generate a file called `giant_height_summary_statistics.tsv.gz` with the same columns as the SNP table and two additional columns:

* `beta_hat`: the estimated effect size via OLS
* `se`: the standard error of the estimated effect size

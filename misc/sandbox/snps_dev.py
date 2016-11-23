import sys
sys.path.append('../../')
import gwass

vcf_path = '/home/jhmarcus/novembre_lab/data/external_public/1kg_phase3/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.biallelic_snps.vcf.gz' 
ref_genome_path = '/home/jhmarcus/novembre_lab/data/external_public/reference_genomes/hs37d5.fa'
snps = gwass.Snps(vcf_path, ref_genome_path)
snps.write_snps()

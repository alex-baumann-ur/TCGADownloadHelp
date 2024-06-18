import pandas as pd

vcf = pd.read_table(snakemake.input[0], comment='#', names=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','TUMOR'])
bed = pd.read_table(snakemake.input[1], compression='gzip')

out_txt = snakemake.output[0]

# VCF


# BED
bed['svclass'].value_counts().to_csv(out_txt, sep='\t')




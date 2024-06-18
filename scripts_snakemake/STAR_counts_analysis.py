import pandas as pd

star = pd.read_table(snakemake.input[0], skiprows=[0,2,3,4,5])

out_file = snakemake.output[0]

star_pseudogene = star[star['gene_type']=='transcribed_unitary_pseudogene']

star_pseudogene.to_csv(out_file, index=False, sep='\t')
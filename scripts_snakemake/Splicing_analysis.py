import pandas as pd

splice = pd.read_table(snakemake.input[0], compression='gzip')

out_file = snakemake.output[0]

splice_chr1 = splice[(splice['#chromosome']=='chr1')&(splice['n_multi_map']>10)]

splice_chr1.to_csv(out_file, index=False, sep='\t')
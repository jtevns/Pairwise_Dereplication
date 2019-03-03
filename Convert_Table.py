import sys 
import pandas as pds

ani_short = pds.read_csv(sys.argv[1],sep="\t")
cov_short = pds.read_csv(sys.argv[2],sep="\t")

ani_long = pds.melt(ani_short,id_vars=['Unnamed: 0'])
cov_long = pds.melt(cov_short,['Unnamed: 0'])

merged_table = pds.merge(ani_long, cov_long,  how='outer', left_on=['Unnamed: 0','variable'],
                            right_on = ['Unnamed: 0','variable'])

merged_table.columns = ['bin1','bin2','ANI','Cov']

merged_table.to_csv("pairwise_long.tsv",sep="\t")

import os
import scipy.stats
import pandas as pd
import sys
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

cors = pd.read_table(sys.argv[1],sep='\t')
cors_sub = pd.DataFrame()
for i in range(len(cors)):
	if(cors.loc[i,'N']>=10):
		str1 = cors.loc[i,'subject']
		str2 = cors.loc[i,'object']
		if(("musicc" in str1) and ("musicc" not in str2)):
			cors_sub = cors_sub.append(cors.iloc[i,:],ignore_index=True)
		elif(("musicc" in str2) and ("musicc" not in str1)):
			cors_sub = cors_sub.append(cors.iloc[i,:],ignore_index=True)
cors_sub.sort_values(by = ["rho","Bonferroni_pval",  "N"], ascending = [False, True,False],ignore_index=True)
cors_sub.to_csv('microbiome_cors_with_analytes', sep='\t', index=False)



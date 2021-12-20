
##To find concept-prefix semantics

import os
import scipy.stats
import pandas as pd
import sys

'''
edges = pd.read_table('wellness_kg_edges_v1.5.tsv', sep='\t')
nodes = pd.read_table('wellness_kg_nodes_v1.5.tsv', sep='\t')
dict_syntax = dict(zip(nodes.loc[:,'id'],nodes.loc[:,'category']))
#for_col = pd.DataFrame({"subject_prefix": [], "subject_semantics": [], "predicate": [], "object_prefix": [], "object_semantics": []})
print("subject_prefix"+'\t'+'subject_semantics'+'\t'+'predicate'+'\t'+'object_prefix'+'\t'+'object_semantics')

for i in range(len(edges)):
	str1 = edges.loc[i,'subject']
	str1 = str1.rpartition(':')[0]
	print(str1,end='\t')
	print(dict_syntax.get(edges.loc[i,'subject']),end='\t')
	print(edges.loc[i,'predicate'],end='\t')
	str2 = edges.loc[i,'object']
	str2 = str2.rpartition(':')[0]
	print(str2,end='\t')
	print(dict_syntax.get(edges.loc[i,'object']))

'''
out = pd.read_table('out', sep='\t')
out = out.drop_duplicates(ignore_index=True)
print("subject_prefix"+'\t'+'subject_semantics'+'\t'+'predicate'+'\t'+'object_prefix'+'\t'+'object_semantics')

for i in range(len(out)):
	print(out.iloc[i,0]+'\t'+out.iloc[i,1]+'\t'+out.iloc[i,2]+'\t'+out.iloc[i,3]+'\t'+out.iloc[i,4])



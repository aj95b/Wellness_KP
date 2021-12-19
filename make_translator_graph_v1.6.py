import os
import scipy.stats
import pandas as pd
import sys

def make_edges():
	current_kg = pd.read_csv('wellness_kg_edges_v1.5.tsv', sep='\t')
	current_kg = current_kg.drop(['weighted_pvalue'],axis=1)
	all_cors = pd.read_csv('join4tables/all_cors.tsv', sep='\t')
	df = pd.read_csv('curies/nodes_new_curies_2021-11-02.tsv', sep='\t')
	cur = dict(zip(df.iloc[:,0],df.iloc[:,1]))

	cols = len(current_kg.columns)
	for i in range(cols):
		if(i==cols-1):
			print(current_kg.columns.values[i])
			break
		print(current_kg.columns.values[i], end='\t')
	

	length = len(current_kg)
	for i in range(length):
		if(current_kg.loc[i,'predicate'] == "biolink:related_to" or not(pd.isna(current_kg.loc[i,'qualifiers']))):
			#vals = current_kg.to_string(index=False, header=False)
			for j in range(cols):
				if(j==cols-1):
					print(current_kg.iloc[i,j])
					break
				print(current_kg.iloc[i,j],end='\t')

	len_new_cors = len(all_cors)
	cols_new = all_cors.columns
	for i in range(len_new_cors):
		if(not(pd.isna(cur.get(all_cors.loc[i,'subject']))) and not(pd.isna(cur.get(all_cors.loc[i,'object'])))):
			print(cur.get(all_cors.loc[i,'subject'])+'\t'+all_cors.loc[i,'predicate']+'\t'+cur.get(all_cors.loc[i,'object'])+'\t'+"RO:0002610"+'\t'+all_cors.loc[i,'subject']+'\t'+all_cors.loc[i,'object']+'\t'+'biolink:Association'+'\t'+str(all_cors.loc[i,'N'])+'\t'+'Spearman Correlation'+'\t'+str(all_cors.loc[i,'rho'])+'\t'+''+'\t'+''+'\t'+str(all_cors.loc[i,'Bonferroni_pval']))


def make_nodes():
	df = pd.read_csv('curies/nodes_new_curies_2021-11-02.tsv', sep='\t')
	cur = dict(zip(df.iloc[:,0],df.iloc[:,1]))
	prote_meta = pd.read_csv('snapshots/arivale_snapshot_ISB_2020-03-16_2156/proteomics_metadata.tsv', sep='\t', comment="#")
	prote_meta_dict = dict(zip(prote_meta.loc[:,'name'], prote_meta.loc[:,'gene_description']))
	print('id'+'\t'+'name'+'\t'+'category'+'\t'+'gene description for UNIPROT ids')
	edges = pd.read_csv('wellness_kg_edges_v1.6.tsv', sep='\t')
	n = len(edges)
	a = []
	for i in range(n):
		a.append(edges.loc[i,'subject_name'])
		a.append(edges.loc[i,'object_name'])
	#print(len(a))
	a = list(set(a))
	#print(len(a))
	#print(cur.get(a[0]))
	
	for i in range(len(a)):
		str1 = cur.get(a[i])
		if(not(pd.isna(str1))):
			print(str1+'\t'+a[i],end='\t')
			if("LOINC" in str1 or "NCIT" in str1):
				print("biolink:ClinicalFinding")
			elif("KEGG.COMPOUND" in str1 or "KEGG.DRUG" in str1 or "HMDB" in str1 or "PUBCHEM" in str1 or "CAS" in str1):
				print("biolink:SmallMolecule")
			elif("UniProt" in str1):
				print("biolink:Protein"+'\t'+prote_meta_dict.get(a[i]))
			elif("KEGG.ORTHOLOGY" in str1):
				print("biolink:MolecularActivity")
			else:
				print("biolink:MolecularEntity")
		
	
		


def main():
	#make_edges()	
	make_nodes()

if __name__=="__main__":
    main()
		


import os
import scipy.stats
import pandas as pd
import sys

def clean_chem(path):
	chem = pd.read_csv(path+'chemistries.tsv', sep='\t', comment="#", dtype={'public_client_id':str})
	#Removing reflexive chemistries
	chem = chem[chem.reflexive != True]
	return chem

def clean_metabolomics(path):
	#Uncomment next line and comment out the line after for snapshots with metabolomics metadata
	met_meta = pd.read_csv(path+'metabolomics_metadata.tsv', sep='\t', comment="#")
	#met_meta = pd.read_csv('snapshots/arivale_snapshot_ISB_2019-05-31_2326/metabolomics_metadata.tsv', sep='\t', comment="#")
	met_meta.columns = met_meta.columns.str.upper()
	#df = pd.read_csv('nodes_new_curies.tsv', sep='\t')
	#nodes_curies_dict = dict(zip(df.iloc[:,0],df.iloc[:,1]))
	v = 'metabolomics.tsv'
	w = 'metabolomics_data.tsv'
	x = 'metabolomics_data_corrected.tsv'
	y = 'metabolomics_corrected.tsv'
	path_met = ''
	path_temp = path
	if os.path.isfile(os.path.join(path_temp,v)):
		path_met = path+v
	elif os.path.isfile(os.path.join(path_temp,w)):
		path_met = path+w
         #snaps.append(path+j)
	elif os.path.isfile(os.path.join(path_temp,x)):
		path_met = path+x
         #snaps.append(path+k)
	elif os.path.isfile(os.path.join(path_temp,y)):
		path_met = path+y
	
	met = pd.read_csv(path_met, sep='\t', comment="#", dtype={'public_client_id':str})
	#Fixing column names metabolites
	if('35' not in met.columns):
		list1 = []
		for k in range(len(met.columns)):
			str1 = met.columns[k]
			str2 = str1.rpartition('.')[-1]
			str2 = str2.replace(':scaled','')
			list1.append(str2)
		new_columns = dict(zip(met.columns, list1))
		met = met.rename(columns=new_columns)
	else:
		new_columns = dict(zip(met_meta['CHEMICAL_ID'].astype(str), met_meta['BIOCHEMICAL_NAME']))
		met = met.rename(columns=new_columns)
	return met

def clean_proteomics(path):
	a = 'proteomics_data.tsv'
	b = 'proteomics_data_corrected.tsv'
	c = 'proteomics_corrected.tsv'
	path_prot = ''
	path_temp = path
	if os.path.isfile(os.path.join(path_temp,a)):
		path_prot = path+a
         #snaps.append(path+j)
	elif os.path.isfile(os.path.join(path_temp,b)):
		path_prot = path+b
         #snaps.append(path+k)
	elif os.path.isfile(os.path.join(path_temp,c)):
		path_prot = path+c

	prot = pd.read_csv(path_prot, sep='\t', comment="#", dtype={'public_client_id':str})
	#Collecting proteins only from healthy people
	prot = sectionalize_prot(prot)
	return prot

def clean_microbiome(path):
	micro = pd.read_csv(path+'microbiome_kos.tsv', sep='\t', comment="#", dtype={'public_client_id':str})
	micro_meta = pd.read_csv(path+'microbiome_analyte_details.tsv', sep='\t', comment="#")
	new_columns = dict(zip(micro_meta['name'].astype(str), micro_meta['display_name_1']))
	micro = micro.rename(columns=new_columns)
	micro = micro.loc[:,~micro.columns.duplicated()]
	return micro

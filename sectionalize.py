import os
import scipy.stats
import pandas as pd
import sys
from clean_data import *


def sectionalize_prot(input_prot):
	prot_healthy = pd.DataFrame()
	drop_prot = []
	for i in range(len(input_prot)):
		for j in range(len(input_prot.columns)):
			if(('CAM' in input_prot.columns[j]) and pd.notnull(input_prot.iloc[i,j])):  #- Cardiometabolic
				drop_prot.append(input_prot.columns[j])
				break
			elif(('CRE' in input_prot.columns[j]) and pd.notnull(input_prot.iloc[i,j])):  #- Cell Regulation
				drop_prot.append(input_prot.columns[j])
				break
			elif(('DEV' in input_prot.columns[j]) and pd.notnull(input_prot.iloc[i,j])):  #- Development
				drop_prot.append(input_prot.columns[j])
				break
			elif(('IMO' in input_prot.columns[j]) and pd.notnull(input_prot.iloc[i,j])):   #- Immuno-Oncology
				drop_prot.append(input_prot.columns[j])
				break
			elif(('IRE'  in input_prot.columns[j]) and pd.notnull(input_prot.iloc[i,j])):   #- Immune Response
				drop_prot.append(input_prot.columns[j])
				break
			elif(('MET' in input_prot.columns[j]) and pd.notnull(input_prot.iloc[i,j])):    #- Metabolism
				drop_prot.append(input_prot.columns[j])
				break
			elif(('NEU1'  in input_prot.columns[j]) and pd.notnull(input_prot.iloc[i,j])):    #- Neurology
				drop_prot.append(input_prot.columns[j])
				break
			elif(('NEX' in input_prot.columns[j]) and pd.notnull(input_prot.iloc[i,j])):        #- Neuro-exploratory
				drop_prot.append(input_prot.columns[j])
				break
			elif(('ODA' in input_prot.columns[j]) and pd.notnull(input_prot.iloc[i,j])):         #- Organ Damage
				drop_prot.append(input_prot.columns[j])
				break
			elif(('ONC2' in input_prot.columns[j]) and pd.notnull(input_prot.iloc[i,j])):        #- Oncology 2
				drop_prot.append(input_prot.columns[j])
				break
			elif(('ONC3' in input_prot.columns[j]) and pd.notnull(input_prot.iloc[i,j]) ):        #- Oncology 3
				drop_prot.append(input_prot.columns[j])
				break
      #'CVD2' #- Cardiovascular 2
      #'CVD3' #- Cardiovascular 3
      #'INF' #- Inflammation
    
			if(j==len(input_prot.columns)-1):
    #print(i)
				prot_healthy = prot_healthy.append(input_prot.iloc[i,:], ignore_index = True)
	for col in (set(drop_prot)):
		prot_healthy = prot_healthy.drop([col], axis=1)
	return prot_healthy




def gender(path_chem, path_prot, path_met, met_meta, nodes_curies_dict, path_clients):

	chem = pd.read_csv(path_chem, sep='\t', comment="#", dtype={'public_client_id':str})
  chem = clean_chem(chem)
	prot = pd.read_csv(path_prot, sep='\t', comment="#", dtype={'public_client_id':str})
  prot = clean_proteomics(prot)
	met = pd.read_csv(path_met, sep='\t', comment="#", dtype={'public_client_id':str})
  met = clean_metabolomics(met)
	clients = pd.read_csv(path_clients, sep='\t', comment="#", dtype={'public_client_id':str})
 


	b = met.merge(prot, how = 'inner', on = ['public_client_id', 'days_in_program'])
	big_table = chem.merge(b, how='inner', on=['public_client_id','days_in_program'])
	big_table = clients.merge(big_table, how='inner', on=['public_client_id'])

	analytes_female = pd.DataFrame()
	analytes_male = pd.DataFrame()
	for i in range(len(big_table)):
		if(big_table.loc[i,'sex'] == 'M'):
			#print(big_table.iloc[i,:])
			analytes_male = analytes_male.append(big_table.iloc[i,:], ignore_index = True)
		elif(big_table.loc[i,'sex'] == 'F'):
			#print(big_table.iloc[i,:])
			analytes_female = analytes_female.append(big_table.iloc[i,:], ignore_index = True)
      
  return analytes_female, analytes_male


def age(path_chem, path_prot, path_met, met_meta, nodes_curies_dict, path_clients):
	chem = pd.read_csv(path_chem, sep='\t', comment="#", dtype={'public_client_id':str})
  chem = clean_chem(chem)
	prot = pd.read_csv(path_prot, sep='\t', comment="#", dtype={'public_client_id':str})
  prot = clean_proteomics(prot)
	met = pd.read_csv(path_met, sep='\t', comment="#", dtype={'public_client_id':str})
  met = clean_metabolomics(met)
	clients = pd.read_csv(path_clients, sep='\t', comment="#", dtype={'public_client_id':str})
 
	b = met.merge(prot, how = 'inner', on = ['public_client_id', 'days_in_program'])
	big_table = chem.merge(b, how='inner', on=['public_client_id','days_in_program'])
	big_table = clients.merge(big_table, how='inner', on=['public_client_id'])

	analytes_young = pd.DataFrame()
	analytes_mid_age = pd.DataFrame()
	analytes_old = pd.DataFrame()
	for i in range(len(big_table)):
		if(big_table.loc[i,'age'] < 35):
			#print(big_table.iloc[i,:])
			analytes_young = analytes_young.append(big_table.iloc[i,:], ignore_index = True)
		elif(big_table.loc[i,'age'] >= 35 and big_table.loc[i,'age'] <= 55):
			#print(big_table.iloc[i,:])
			analytes_mid_age = analytes_mid_age.append(big_table.iloc[i,:], ignore_index = True)
		elif(big_table.loc[i,'age'] > 55):
			analytes_old = analytes_old.append(big_table.iloc[i,:], ignore_index = True)

  return analytes_young, analytes_mid_age, analytes_old



def ancestry(path_chem, path_prot, path_met, met_meta, nodes_curies_dict, path_clients):
        chem = pd.read_csv(path_chem, sep='\t', comment="#", dtype={'public_client_id':str})
        chem = clean_chem(chem)
	      prot = pd.read_csv(path_prot, sep='\t', comment="#", dtype={'public_client_id':str})
        prot = clean_proteomics(prot)
	      met = pd.read_csv(path_met, sep='\t', comment="#", dtype={'public_client_id':str})
        met = clean_metabolomics(met)
	      clients = pd.read_csv(path_clients, sep='\t', comment="#", dtype={'public_client_id':str})

        b = met.merge(prot, how = 'inner', on = ['public_client_id', 'days_in_program'])
        big_table = chem.merge(b, how='inner', on=['public_client_id','days_in_program'])
        big_table = clients.merge(big_table, how='inner', on=['public_client_id'])
	
        analytes_white = pd.DataFrame()
        analytes_black = pd.DataFrame()
        analytes_spanish = pd.DataFrame()
        analytes_east_asian = pd.DataFrame()
        analytes_south_asian = pd.DataFrame()
        analytes_other = pd.DataFrame()
        analytes_midd_eastern = pd.DataFrame()
        analytes_native_american = pd.DataFrame()
	
        for i in range(len(big_table)):
                if(big_table.loc[i,'race'] == 'white' or big_table.loc[i,'race'] == 'ashkenazi jewish'):
                        #print(big_table.iloc[i,:])
                        analytes_white = analytes_white.append(big_table.iloc[i,:], ignore_index = True)
                elif(big_table.loc[i,'race'] == 'black or african-american' or big_table.loc[i,'race'] == 'afro-caribbean'):
			#print(big_table.iloc[i,:])
                        analytes_black = analytes_black.append(big_table.iloc[i,:], ignore_index = True)
                elif(big_table.loc[i,'race'] == 'hispanic latino or spanish origin'):
                        analytes_spanish = analytes_spanish.append(big_table.iloc[i,:], ignore_index = True)
                elif(big_table.loc[i,'race'] == 'east asian'):
                        analytes_east_asian = analytes_east_asian.append(big_table.iloc[i,:], ignore_index = True)
                elif(big_table.loc[i,'race'] == 'south asian'):
                        analytes_south_asian = analytes_south_asian.append(big_table.iloc[i,:], ignore_index = True)
                elif(big_table.loc[i,'race'] == 'other'):
                        analytes_other = analytes_other.append(big_table.iloc[i,:], ignore_index = True)
                elif(big_table.loc[i,'race'] == 'middle eastern or north african' or big_table.loc[i,'race'] == 'sephardic jewish'):
                        analytes_midd_eastern = analytes_midd_eastern.append(big_table.iloc[i,:], ignore_index = True)
                elif(big_table.loc[i,'race'] == 'native hawaiian or other pacific islander' or 'american indian or alaska native'):
                        analytes_native_american = analytes_native_american.append(big_table.iloc[i,:], ignore_index = True)
                    
        return analytes_white, analytes_black, analytes_spanish, analytes_east_asian, analytes_south_asian, analytes_other, analytes_midd_eastern, analytes_native_american

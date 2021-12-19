#Computes regressors of analytes that have maximum significant correlations with other analytes
#using ridge regression: primary reason- no. of variables >> observations
#Parameter 1: Path to data snapshot
#Parameter 2: Path to latest correlations
import pandas as pd
import numpy as np
from math import log
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from regressors import stats

#model for predicting IL-6, only contains analytes that are correlated with IL-6
def for_IL6_cors(big_table, well_edges):
	big_table_2 = big_table.copy()
	
	analytes_of_IL6 = []
	for i in range(len(well_edges)):
		if(well_edges.loc[i,'subject']=='UniProtKB:P05231' or well_edges.loc[i,'subject']=='LOINC:26881-3'):
			analytes_of_IL6.append(well_edges.loc[i,'object_name'])
		elif(well_edges.loc[i,'object']=='UniProtKB:P05231' or well_edges.loc[i,'object']=='LOINC:26881-3'):
			analytes_of_IL6.append(well_edges.loc[i,'subject_name'])


	IL6 = ['IL-6']
	unionIL6 = sorted(list(set(analytes_of_IL6)) + list(set(IL6)))
	data_IL6 = pd.DataFrame(columns = unionIL6)

	for i in range(len(big_table_2)):
		for col in (set(unionIL6) & set(big_table_2.columns)):
			data_IL6.loc[i,col] = big_table_2.loc[i,col]

	y = data_IL6.loc[0:,'CVD2_P05231']
	y = y.fillna(y.mean())
	y = np.array(y, dtype=float)
	data_IL6.drop(['INF_P05231','CVD2_P05231', 'IRE_P05231'], inplace=True, axis=1)
	data_IL6 = data_IL6.fillna(data_IL6.mean())
	
	clf1 = Ridge(alpha=1.0)
	clf1.fit(data_IL6, y)
	Ridge(normalize=True)
	clf1.coef_
	a = np.array(clf1.coef_,dtype=float)
	b_pvals = stats.coef_pval(clf1, data_IL6, y)
	b = pd.DataFrame({"analytes": [], "ridge_coef": [], "p_vals": []})
	for i in range(len(a)):
		b.loc[i,'ridge_coef'] = a[i]
		b.loc[i,'analytes'] = big_table_2.columns[i]
		b.loc[i, 'p_vals'] = b_pvals[i]
	b = b.sort_values(by=['ridge_coef'],ignore_index=True)
	return b
	#b.to_csv('cors_IL6_coef.tsv',sep='\t', index=False)

#model for IL-6, containing all other analytes 
def for_IL6_all_data(big_table, well_edges):
	big_table_1 = big_table.copy()
	y = big_table_1.loc[0:,'CVD2_P05231']
	y = y.fillna(y.mean())
	y = np.array(y, dtype=float)
	big_table_1.drop(['INF_P05231','CVD2_P05231', 'IRE_P05231', 'lidocaine', 'X - 21315','cotinine N-oxide'], inplace=True, axis=1)
	big_table_1 = big_table_1.fillna(big_table_1.mean())
	clf = Ridge(alpha=1.0)
	clf.fit(big_table_1, y)
	Ridge(normalize=True)
	clf.coef_
	a = np.array(clf.coef_,dtype=float)
	b_pvals = stats.coef_pval(clf, big_table_1, y)
	b = pd.DataFrame({"analytes": [], "ridge_coef": [], "p_vals": []})
	for i in range(len(a)):
		b.loc[i,'ridge_coef'] = a[i]
		b.loc[i,'analytes'] = big_table_1.columns[i]
		b.loc[i, 'p_vals'] = b_pvals[i]
	b = b.sort_values(by=['ridge_coef'],ignore_index=True)
	#return b
	#b.to_csv('all_data_regressors_IL6_coef.tsv',sep='\t', index=False)

#model for Adiponectin, containing all other analytes
def for_Adipo_serum_all_data(big_table, well_edges):
	big_table_1 = big_table.copy()
	y = big_table_1.loc[0:,'ADIPONECTIN, SERUM']
	y = y.fillna(y.mean())
	y = np.array(y, dtype=float)
	big_table_1.drop(['ADIPONECTIN, SERUM', 'lidocaine', 'X - 21315','cotinine N-oxide'], inplace=True, axis=1)
	big_table_1 = big_table_1.fillna(big_table_1.mean())
	clf = Ridge(alpha=1.0)
	clf.fit(big_table_1, y)
	Ridge(normalize=True)
	clf.coef_
	a = np.array(clf.coef_,dtype=float)
	b_pvals = stats.coef_pval(clf, big_table_1, y)
	b = pd.DataFrame({"analytes": [], "ridge_coef": [], "p_vals": []})
	for i in range(len(a)):
		b.loc[i,'ridge_coef'] = a[i]
		b.loc[i,'analytes'] = big_table_1.columns[i]
		b.loc[i, 'p_vals'] = b_pvals[i]
	b = b.sort_values(by=['ridge_coef'],ignore_index=True)
	#return b
	#b.to_csv('all_data_regressors_Adipo_serum_coef.tsv',sep='\t', index=False)


def for_Adipo_serum_cors(big_table, well_edges):
	big_table_2 = big_table.copy()
	
	analytes_of_adipo = []
	for i in range(len(well_edges)):
		if(well_edges.loc[i,'subject']=='LOINC:47828-9'):
			analytes_of_adipo.append(well_edges.loc[i,'object_name'])
		elif(well_edges.loc[i,'object']=='LOINC:47828-9'):
			analytes_of_adipo.append(well_edges.loc[i,'subject_name'])

	union = sorted(list(set(analytes_of_adipo)))
	data_adipo = pd.DataFrame(columns = union)

	for i in range(len(big_table_2)):
		for col in (set(union) & set(big_table_2.columns)):
			data_adipo.loc[i,col] = big_table_2.loc[i,col]

	y = big_table_2.loc[0:,'ADIPONECTIN, SERUM']
	y = y.fillna(y.mean())
	y = np.array(y, dtype=float)
	#data_adipo.drop(['ADIPONECTIN, SERUM'], inplace=True, axis=1)
	data_adipo = data_adipo.fillna(data_adipo.mean())

	clf1 = Ridge(alpha=1.0)
	clf1.fit(data_adipo, y)
	Ridge(normalize=True)
	clf1.coef_
	a = np.array(clf1.coef_,dtype=float)
	b_pvals = stats.coef_pval(clf1, data_adipo, y)
	b = pd.DataFrame({"analytes": [], "ridge_coef": [], "p_vals": []})
	for i in range(len(a)):
		b.loc[i,'ridge_coef'] = a[i]
		b.loc[i,'analytes'] = big_table_2.columns[i]
		b.loc[i, 'p_vals'] = b_pvals[i]
	b = b.sort_values(by=['ridge_coef'],ignore_index=True)
	return b

def for_arachidoyl_GPC_C04230_all_data(big_table, well_edges): # a metabolite
	big_table_1 = big_table.copy()
	y = big_table_1.loc[0:,'1-arachidoyl-GPC (20:0)']
	y = y.fillna(y.mean())
	y = np.array(y, dtype=float)
	big_table_1.drop(['1-arachidoyl-GPC (20:0)', 'lidocaine', 'X - 21315','cotinine N-oxide'], inplace=True, axis=1)
	big_table_1 = big_table_1.fillna(big_table_1.mean())
	clf = Ridge(alpha=1.0)
	clf.fit(big_table_1, y)
	Ridge(normalize=True)
	clf.coef_
	a = np.array(clf.coef_,dtype=float)
	b_pvals = stats.coef_pval(clf, big_table_1, y)
	b = pd.DataFrame({"analytes": [], "ridge_coef": [], "p_vals": []})
	for i in range(len(a)):
		b.loc[i,'ridge_coef'] = a[i]
		b.loc[i,'analytes'] = big_table_1.columns[i]
		b.loc[i, 'p_vals'] = b_pvals[i]
	b = b.sort_values(by=['ridge_coef'],ignore_index=True)
	#return b
	#b.to_csv('all_data_regressors_arachidoyl_GPC_C04230_coef.tsv',sep='\t', index=False)

def for_arachidoyl_GPC_C04230_cors(big_table, well_edges):
	big_table_2 = big_table.copy()
	
	analytes_of = []
	for i in range(len(well_edges)):
		if(well_edges.loc[i,'subject_name']=='1-arachidoyl-GPC (20:0)'):
			analytes_of.append(well_edges.loc[i,'object_name'])
		elif(well_edges.loc[i,'object_name']=='1-arachidoyl-GPC (20:0)'):
			analytes_of.append(well_edges.loc[i,'subject_name'])

	union = sorted(list(set(analytes_of)))
	data = pd.DataFrame(columns = union)

	for i in range(len(big_table_2)):
		for col in (set(union) & set(big_table_2.columns)):
			data.loc[i,col] = big_table_2.loc[i,col]

	y = big_table_2.loc[0:,'1-arachidoyl-GPC (20:0)']
	y = y.fillna(y.mean())
	y = np.array(y, dtype=float)
	#data.drop(['1-arachidoyl-GPC (20:0)', '1-myristoyl-GPC (14:0)', '1-margaroyl-GPC (17:0)'], inplace=True, axis=1)
	data = data.fillna(data.mean())

	clf1 = Ridge(alpha=1.0)
	clf1.fit(data, y)
	Ridge(normalize=True)
	clf1.coef_
	a = np.array(clf1.coef_,dtype=float)
	b_pvals = stats.coef_pval(clf1, data, y)
	b = pd.DataFrame({"analytes": [], "ridge_coef": [], "p_vals": []})
	for i in range(len(a)):
		b.loc[i,'ridge_coef'] = a[i]
		b.loc[i,'analytes'] = big_table_2.columns[i]
		b.loc[i, 'p_vals'] = b_pvals[i]
	b = b.sort_values(by=['ridge_coef'],ignore_index=True)
	return b

def for_glomerular_fr_cors(big_table, well_edges): #a clinical lab result
	big_table_2 = big_table.copy()
	
	analytes_of = []
	for i in range(len(well_edges)):
		if(well_edges.loc[i,'subject']=='LOINC:62238-1'):
			analytes_of.append(well_edges.loc[i,'object_name'])
		elif(well_edges.loc[i,'object']=='LOINC:62238-1'):
			analytes_of.append(well_edges.loc[i,'subject_name'])


	
	union = sorted(list(set(analytes_of)))
	data = pd.DataFrame(columns = union)

	for i in range(len(big_table_2)):
		for col in (set(union) & set(big_table_2.columns)):
			data.loc[i,col] = big_table_2.loc[i,col]

	y = big_table_2.loc[0:,'GFR, MDRD']
	y = y.fillna(y.mean())
	y = np.array(y, dtype=float)
	data.drop(['GFR, MDRD, AFRICAN AM'], inplace=True, axis=1)
	data = data.fillna(data.mean())

	clf1 = Ridge(alpha=1.0)
	clf1.fit(data, y)
	Ridge(normalize=True)
	clf1.coef_
	a = np.array(clf1.coef_,dtype=float)
	b_pvals = stats.coef_pval(clf1, data, y)
	b = pd.DataFrame({"analytes": [], "ridge_coef": [], "p_vals": []})
	for i in range(len(a)):
		b.loc[i,'ridge_coef'] = a[i]
		b.loc[i,'analytes'] = big_table_2.columns[i]
		b.loc[i, 'p_vals'] = b_pvals[i]
	b = b.sort_values(by=['ridge_coef'],ignore_index=True)
	return b
	#b.to_csv('cors_IL6_coef.tsv',sep='\t', index=False)


def for_glomerular_fr_all_data(big_table, well_edges):
	big_table_1 = big_table.copy()
	y = big_table_1.loc[0:,'GFR, MDRD']
	y = y.fillna(y.mean())
	y = np.array(y, dtype=float)
	big_table_1.drop(['GFR, MDRD','GFR, MDRD, AFRICAN AM', 'lidocaine', 'X - 21315','cotinine N-oxide'], inplace=True, axis=1)
	big_table_1 = big_table_1.fillna(big_table_1.mean())
	clf = Ridge(alpha=1.0)
	clf.fit(big_table_1, y)
	Ridge(normalize=True)
	clf.coef_
	a = np.array(clf.coef_,dtype=float)
	b_pvals = stats.coef_pval(clf, big_table_1, y)
	b = pd.DataFrame({"analytes": [], "ridge_coef": [], "p_vals": []})
	for i in range(len(a)):
		b.loc[i,'ridge_coef'] = a[i]
		b.loc[i,'analytes'] = big_table_1.columns[i]
		b.loc[i, 'p_vals'] = b_pvals[i]
	b = b.sort_values(by=['ridge_coef'],ignore_index=True)
	return b
	#b.to_csv('all_data_regressors_glomerular_fr_coef.tsv',sep='\t', index=False)


def for_crp_all_data(big_table, well_edges):#CRP: C-Reactive protein
	big_table_1 = big_table.copy()
	y = big_table_1.loc[0:,'CRP HIGH SENSITIVITY']
	y = y.fillna(y.mean())
	y = np.array(y, dtype=float)
	big_table_1.drop(['CRP HIGH SENSITIVITY', 'lidocaine', 'X - 21315','cotinine N-oxide'], inplace=True, axis=1)
	big_table_1 = big_table_1.fillna(big_table_1.mean())
	clf = Ridge(alpha=1.0)
	clf.fit(big_table_1, y)
	Ridge(normalize=True)
	clf.coef_
	a = np.array(clf.coef_,dtype=float)
	b_pvals = stats.coef_pval(clf, big_table_1, y)
	b = pd.DataFrame({"analytes": [], "ridge_coef": [], "p_vals": []})
	for i in range(len(a)):
		b.loc[i,'ridge_coef'] = a[i]
		b.loc[i,'analytes'] = big_table_1.columns[i]
		b.loc[i, 'p_vals'] = b_pvals[i]
	b = b.sort_values(by=['ridge_coef'],ignore_index=True)
	return b
	#b.to_csv('all_data_regressors_crp_coef.tsv',sep='\t', index=False)

def for_crp_cors(big_table, well_edges):
	big_table_2 = big_table.copy()
	
	analytes_of_crp = []
	for i in range(len(well_edges)):
		if(well_edges.loc[i,'subject']=='LOINC:30522-7'):
			analytes_of_crp.append(well_edges.loc[i,'object_name'])
		elif(well_edges.loc[i,'object']=='LOINC:30522-7'):
			analytes_of_crp.append(well_edges.loc[i,'subject_name'])

	union = sorted(list(set(analytes_of_crp)))
	data_crp = pd.DataFrame(columns = union)

	for i in range(len(big_table_2)):
		for col in (set(union) & set(big_table_2.columns)):
			data_crp.loc[i,col] = big_table_2.loc[i,col]

	y = big_table_2.loc[0:,'CRP HIGH SENSITIVITY']
	y = y.fillna(y.mean())
	y = np.array(y, dtype=float)
	#data_adipo.drop(['ADIPONECTIN, SERUM'], inplace=True, axis=1)
	data_crp = data_crp.fillna(data_crp.mean())

	clf1 = Ridge(alpha=1.0)
	clf1.fit(data_crp, y)
	Ridge(normalize=True)
	clf1.coef_
	a = np.array(clf1.coef_,dtype=float)
	b_pvals = stats.coef_pval(clf1, data_crp, y)
	b = pd.DataFrame({"analytes": [], "ridge_coef": [], "p_vals": []})
	for i in range(len(a)):
		b.loc[i,'ridge_coef'] = a[i]
		b.loc[i,'analytes'] = big_table_2.columns[i]
		b.loc[i, 'p_vals'] = b_pvals[i]
	b = b.sort_values(by=['ridge_coef'],ignore_index=True)
	return b


def non_cr(all_rgr, cor_rgr):
	merged = pd.merge(all_rgr, cor_rgr, how='inner', on=['analytes'])
	result =all_rgr[~all_rgr.analytes.isin(merged.analytes)]
	result = result.sort_values(by=['ridge_coef'],ignore_index=True)
	return result

def main():
  path_to_snapshot = sys.argv[1]
  path_to_correlations = sys.argv[2]
  chem = pd.read_csv(path_to_snapshot, sep='\t', comment="#", dtype={'public_client_id':str})
  chem = clean_chem(chem)
  prot = pd.read_csv(path_to_snapshot, sep='\t', comment="#", dtype={'public_client_id':str})
  prot = clean_proteomics(prot)
  met = pd.read_csv(path_to_snapshot, sep='\t', comment="#", dtype={'public_client_id':str})
  met = clean_metabolomics(met)

  

  b = met.merge(prot, how = 'inner', on = ['public_client_id', 'days_in_program'])
  big_table = chem.merge(b, how='inner', on=['public_client_id','days_in_program'])
  well_edges = pd.read_csv(path_to_correlations, sep='\t')

  drop_list = ['sample_id', 'sample_id_x', 'sample_id_y', 'SAMPLE_ID', 'metabolon_collection_date', 'days_since_first_call', 'days_since_first_call_x','days_since_first_call_y',
              'days_since_first_draw', 'days_since_first_draw_x','days_since_first_draw_y','month', 'month_x','month_y','weekday', 'weekday_x','weekday_y','public_client_id','days_in_program','days_in_program_x','days_in_program_y','Chip_ID_CVD2', 'Chip_ID_CVD3', 'Chip_ID_INF','CVD3_Q99727', 'Chip_ID_CAM' ,'Chip_ID_CRE','Chip_ID_DEV' ,'Chip_ID_IRE','Chip_ID_MET' ,'Chip_ID_NEU1', 'Chip_ID_NEX','Chip_ID_ODA','Chip_ID_ONC2','Chip_ID_ONC3','INF_P01732','batch_CVD2', 'batch_CVD3', 'batch_INF', 'verbose_batch_CVD2', 'verbose_batch_CVD3', 'verbose_batch_INF', 'vendor','vendor_x','vendor_y', 'season', 'season_x','season_y', 'vendor','vendor_observation_id','observation_id','reflexive','month','fasting','vobs_id','accession_number','panel_code', 'group','source_id','status','collection_date_BSTG', 'mrn_BSTG','public_client_id_BSTG', 'vobs_id_BSTG', 'genome_id_BSTG', 'accession_number_BSTG', 'how_joined' ,'merge_priority']

  big_table.set_index('public_client_id', 'days_in_program')
  big_table = big_table.groupby(['public_client_id'],as_index=False).first()
  for col in (set(drop_list) & set(big_table.columns)):
        big_table = big_table.drop([col], axis=1)



  #all_regressors =
  for_IL6_all_data(big_table, well_edges)
  #correlated_regressors = for_IL6_cors(big_table, well_edges)
  #all_regressors = 
  for_Adipo_serum_all_data(big_table, well_edges)
  #correlated_regressors = for_Adipo_serum_cors(big_table, well_edges)
  #all_regressors = 
  for_arachidoyl_GPC_C04230_all_data(big_table, well_edges)
  #correlated_regressors = for_arachidoyl_GPC_C04230_cors(big_table, well_edges)
  #all_regressors = 
  for_glomerular_fr_all_data(big_table, well_edges)
  #correlated_regressors = for_glomerular_fr_cors(big_table, well_edges)
  #all_regressors = 
  for_crp_all_data(big_table, well_edges)
  #correlated_regressors = for_crp_cors(big_table, well_edges)
  #non_correlated_regressors = non_cr(all_regressors, correlated_regressors)
  #non_correlated_regressors.to_csv('non_correlated_regressors_crp.tsv',sep='\t', index=False)
  
if __name__=="__main__":
    main()

#@Author: Arpita Joshi. Acknowledging help from Andrew Magis, ISB

#Parameter 1: path to snapshot
#Parameter 2: path to latest correlations of the corresponding snapshot
#Parameters 3 and 4: the two tables whose joint correlation is received via Parameter #2
#Parameter 5: type of merge (inner, outer, right or left), correlations of variables from one table might already exist because of joint correlation with
            # a different table. Would save compute time to not use outer (the union) each time

#######Notes: to extract coefficients and p-values as they get written in a way that has both tabs and spaces in the output##########
#### x_coef = x.loc[:,'intercept_coef'] x_pval = x.loc[:,'intervention_coef']
#### x_pval = x_pval.to_frame()
#### x_pval = pd.DataFrame(x_pval.intervention_coef.str.split(' ',3).tolist(),columns=['a','b','c','d'])
#### x_pval = x_pval.astype(float)
#### x_pval.sort_values('d')

import os
import scipy.stats
import pandas as pd
import sys
from cors4table import *
from sectionalize import sectionalize_prot
pd.set_option('display.max_columns', None)
import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.sandbox.stats.multicomp
from statsmodels.genmod.families import family, links
import itertools

def compute_interactions(sub):
	col_1 = sub.columns[0]
	col_2 = sub.columns[1]
	col_3 = sub.columns[2]
	sub = sub.rename(columns={col_1:'analyte1', col_2:'analyte2', col_3:'analyte3'})
	family_type = family.Gaussian()
	family_type.link = links.identity()
	family_name = 'Gaussian'
	family_link = 'Identity'
	if (sub.loc[:,'analyte1'].skew() > 1.5) or (sub.loc[:,'analyte1'].skew() < -1.5):

            #logger.info('Setting gamma family for skewed analyte %s'%(col))

            # Set any zero values to 1/2 the smallest value
		sub.loc[sub['analyte1']==0, 'analyte1'] = (sub.loc[sub['analyte1']>0, 'analyte1'].min() / 2.0)

		family_type = family.Gamma()
		family_type.link = links.log()
		family_name = 'Gamma'
		family_link = 'Log'
	
	
	try:
		ols_model = 'analyte1 ~ analyte3*analyte2'
		fitted_model = smf.glm(ols_model, data=sub, family=family_type, missing='drop').fit(maxiter=2000)
		if(fitted_model.pvalues.iloc[-1] <= 0.05):
			print(col_1+'\t'+col_2+'\t'+col_3+'\t'+ str(len(fitted_model.fittedvalues))+'\t'+str(fitted_model.converged), end='\t')#+'\t'+*fitted_model.params+'\t'+ *fitted_model.pvalues)
			print(*fitted_model.params, end='\t')
			print(*fitted_model.pvalues)
			#print(fitted_model.summary())
	except Exception as e:
		#logger.info('Failed analytes {} {} with error {}'.format(col1, col2, str(e)))
		pass
	#print(fitted_model.summary())
	
def main():
	path_to_snapshot = sys.argv[1]
	table_1 = sys.argv[3]
	table_2 = sys.argv[4]
	table_frame_1 = pd.DataFrame()
	table_frame_2 = pd.DataFrame()
	
	if (table_1=='chemistries'):
		table_frame_1 = clean_chem(path_to_snapshot)
	elif(table_1=='metabolomics'):
		table_frame_1 = clean_metabolomics(path_to_snapshot)
	elif(table_1=='proteomics'):
		table_frame_1 = clean_proteomics(path_to_snapshot)
	elif(table_1=='microbiome'):
		table_frame_1 = clean_microbiome(path_to_snapshot)
	
	if (table_2=='chemistries'):
		table_frame_2 = clean_chem(path_to_snapshot)
	elif(table_2=='metabolomics'):
		table_frame_2 = clean_metabolomics(path_to_snapshot)
	elif(table_2=='proteomics'):
		table_frame_2 = clean_proteomics(path_to_snapshot)
	elif(table_2=='microbiome'):
		table_frame_2 = clean_microbiome(path_to_snapshot)
		

	
	type_of_merge = sys.argv[5]
	data_table = pd.merge(table_frame_1, table_frame_2, how=type_of_merge, on=['public_client_id','days_in_program'])
	
	drop_list =['sample_id','sample_id_x','sample_id_y','sample_id_z','SAMPLE_ID','metabolon_collection_date','days_since_first_call','days_since_first_call_x',
'days_since_first_call_y','days_since_first_call_z','days_since_first_draw','days_since_first_draw_x','days_since_first_draw_y','days_since_first_draw_z',
'month','month_x','month_y','month_z','weekday','weekday_x','weekday_y','weekday_z','public_client_id','days_in_program','days_in_program_x','days_in_program_y',
'days_in_program_z','Chip_ID_CVD2','Chip_ID_CVD3','Chip_ID_INF','CVD3_Q99727','Chip_ID_CAM' ,'Chip_ID_CRE','Chip_ID_DEV' ,'Chip_ID_IRE','Chip_ID_MET' ,
'Chip_ID_NEU1','Chip_ID_NEX','Chip_ID_ODA','Chip_ID_ONC2','Chip_ID_ONC3','INF_P01732','batch_CVD2','batch_CVD3','batch_INF','verbose_batch_CVD2',
'verbose_batch_CVD3','verbose_batch_INF','vendor','vendor_x','vendor_y','vendor_z','season','season_x','season_y','season_z','vendor','vendor_x','vendor_y',
'vendor_observation_id','vendor_observation_id_x','vendor_observation_id_y',
'observation_id','reflexive','month','fasting','vobs_id','accession_number','panel_code','group','source_id','status',
'collection_date_BSTG','mrn_BSTG','public_client_id_BSTG','vobs_id_BSTG','genome_id_BSTG','accession_number_BSTG','how_joined' ,'merge_priority','genome_id',
'genome_vendor','multiple_genome_ids','client_id','user_id','race','has_research_consent','enterprise','coach','region','country','current_program',
'latest_program','is_helix','age']
	data_table = table_frame_1
	data_table = data_table.set_index('public_client_id', 'days_in_program')
	data_table = data_table.groupby(['public_client_id'],as_index=False).first()

	for col in (set(drop_list) & set(data_table.columns)):
		data_table = data_table.drop([col], axis=1)

	
	print('subject'+'\t'+'object'+'\t'+'interacting_analyte(IA)'+'\t'+'N'+'\t'+'converged'+'\t'+'intercept_coef'+
'\t'+'IA_coef'+'\t'+'object_coef'+'\t'+'interaction_coef(IA:object)'+'\t'+'intercept_pval'+'\t'+'IA_pval'+'\t'+
'object_pval'+'\t'+'interaction_pval(IA:object)')
	
	for i in range(len(cors)):
		if(cors.loc[i,'rho']>=0.97): #ignore correlations that are likely the same analyte measured twice
			continue

		
		for j in range(len(data_table.columns)):
			#include more stipulations as necessary, to exclude interaction that make no sense
			if(data_table.columns[j]==cors.loc[i,'subject'] or data_table.columns[j]==cors.loc[i,'object'] or data_table.columns[j] not in cors.values):#or data_table.columns[j] not in micro.columns): 
				continue
			
			sub = data_table[[cors.loc[i,'subject'], cors.loc[i,'object'], data_table.columns[j]]].copy()				
			compute_interactions(sub)
	
  
if __name__=="__main__":
    main()

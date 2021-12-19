#@Author: Arpita Joshi
#Parameter 1: path to snapshot
#Parameter 2 and 3: the two tables whose joint correlation is desired
#Parameter 4: type of merge (inner, outer, right or left), correlations of variables from one table might already exist because of joint correlation with a different table. Would save compute time to not use outer (the union) each time.

import os
import scipy.stats
import pandas as pd
import sys

from sectionalize import sectionalize_prot


def calc_cors(table):
	num_analytes = len(table.columns)
	table_size = len(table)
	bon_coeff = num_analytes*(num_analytes-1)/2
	#Calculate only for unique curie IDs / concept pairs
	print("subject"+'\t'+"predicate"+'\t'+"object"+'\t'+"N"+'\t'+"rho"+'\t'+"Bonferroni_pval")
	for i in range(num_analytes):
		for j in range(num_analytes):
                        if(j<=i):
                                continue
			
                        #Evaluate concept pair counts
                        num_not_null = 0
                        for a in range(table_size):
                                if(pd.notnull(table.iloc[a,i]) and pd.notnull(table.iloc[a,j])):
                                        num_not_null+=1

                        if(num_not_null<10):
                                continue

                        corr_pval = scipy.stats.spearmanr(table.iloc[:,i], table.iloc[:,j], nan_policy='omit')
                        if(bon_coeff * corr_pval[1] <= 0.05):
                                print(table.columns[i], end='\t')
                                print("correlated_with", end='\t')
                                print(table.columns[j],end='\t')
                                print(num_not_null,end='\t')
                                print(corr_pval[0],end='\t')
                                print(corr_pval[1])
		

def main():
	path_to_snapshot = sys.argv[1]
	table_1 = sys.argv[2]
	#table_2 = sys.argv[3]
	table_frame_1 = pd.DataFrame()
	#table_frame_2 = pd.DataFrame()
	
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
		
	type_of_merge = sys.argv[4]
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

	
	calc_cors(data_table)
	
  
if __name__=="__main__":
    main()

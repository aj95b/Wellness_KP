#@Author: Arpita Joshi. Acknowledging help from Andrew Magis, ISB

#Parameter 1: path to snapshot
#Parameter 2: path to latest correlations of the corresponding snapshot

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
	#sub = sub.rename(columns={sub.columns[1]:'analyte2'})
	#sub = sub.rename(columns={sub.columns[2]:'analyte3'})
	#print(sub.head())
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
	#ols_model = 'analyte1 ~ age*analyte2*C(sex) + C(season) + BMI_CALC + C(meds_cholesterol) + C(meds_blood_sugar) + C(meds_blood_pressure) + C(vendor) + PC1 + PC2 + PC3 + PC4'
	
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

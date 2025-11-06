#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os 
import sys 
import glob
import gc
import pandas as pd
import numpy as np 
import re

import statsmodels.api as sm
import statsmodels.formula.api as smf


vartype = dict.fromkeys(['r2', 'Sigma', 'R1t2', 'deltat', 'it', 't', 'N_SSS', 'N_SSI', 'N_SSR',
       'N_SIS', 'N_SIR', 'N_SRS', 'N_SRI', 'N_SRR', 'N_ISS', 'N_ISR', 'N_IRS',
       'N_IRR', 'N_RSS', 'N_RSI', 'N_RSR', 'N_RIS', 'N_RIR', 'N_RRS', 'N_RRI',
       'N_RRR', 'Inc_SSI', 'Inc_SIS', 'Inc_SIR', 'Inc_SRI', 'Inc_ISS',
       'Inc_ISR', 'Inc_IRS', 'Inc_IRR', 'Inc_RSI', 'Inc_RIS', 'Inc_RIR',
       'Inc_RRI'],'int32')

vartype['r2'] = 'str'
vartype['Sigma'] = 'str'
vartype['R1t2'] = 'str'


def model(df):
    glm_binom = smf.glm('Inc2 + Inc1 ~ t_from_emerge',df, family=sm.families.Binomial())
    result=glm_binom.fit()
    return result

def ModelSumm(result):
    summ = pd.concat([result.params.rename('coeff'),result.bse.rename('std_err'),result.pvalues.rename('pvalue'),
             result.conf_int().rename({0:'ci_2.5',1:'ci_97.5'},axis=1)],axis=1)
    summ = summ.reset_index().melt(id_vars='index').set_index(['index','variable']).squeeze()

    return summ

def ModelSumm2(result):
    summ = pd.Series({'count':result.nobs,'pseudo_rsquared':result.pseudo_rsquared(),'pearson_chi2':result.pearson_chi2,
                     'llf':result.llf,'llnull':result.llnull,'deviance':result.deviance,
                     'null_deviance':result.null_deviance,'aic':result.aic})
    return summ


path = sys.argv[1]
print('Hi This is the path: ', path, flush=True)
FILES = glob.glob(path+'TimeSerie_*.csv', recursive=True)
FILES = pd.DataFrame(FILES,columns=['Filename'])
FILES['Net'] = FILES['Filename'].str.extract('(_Net\d+)', expand=False).str.replace('_','')
#FILES['Net'] = FILES['Filename'].str.extract('(Net_.._r)', expand=False).str.replace('_r','')
FILES['Sigma2'] = FILES['Filename'].str.extract('(sigma2=\d+\.?\d*)', expand=False).str.replace('sigma2=','')
FILES['Netrealisation'] = FILES['Filename'].str.extract('(r\d+_)', expand=False).str.replace('_','').str.replace('r','')
FILES['Netrealisation'] = FILES['Netrealisation'].astype(int)
FILES.sort_values(by=['Net','Netrealisation'],inplace=True)
FILES.reset_index(inplace=True,drop=True)
# FILES
#print(FILES, flush=True)
Models_All = pd.DataFrame()

Nets = pd.unique(FILES['Net'])

print('Nets are: ', Nets, flush=True)

for net in Nets:
    FILEStemps = FILES[FILES['Net']==net]
    TimeSeDT = pd.DataFrame()
    for i in range(len(FILEStemps)):
        FILE = FILEStemps.iloc[i]['Filename']
        r = FILEStemps.iloc[i]['Netrealisation']
        Sigma2 = FILEStemps.iloc[i]['Sigma2']
        TimeSeDT_temp = pd.read_csv(FILE,dtype = vartype)
        TimeSeDT_temp['Netrealisation'] = r
        TimeSeDT_temp['Sigma2'] = Sigma2
        TimeSeDT_temp['Net'] = net
        TimeSeDT = pd.concat([TimeSeDT, TimeSeDT_temp])
    
    del TimeSeDT_temp
    print('Net:',net, ' ',i,' read tables', flush=True)
    #print(TimeSeDT,flush=True)
    #Models_All = pd.DataFrame()

    TimeSeDT['I1'] = TimeSeDT['N_ISS'] + TimeSeDT['N_ISR'] + TimeSeDT['N_IRS'] + TimeSeDT['N_IRR']
    TimeSeDT['I2'] = TimeSeDT['N_SIS'] + TimeSeDT['N_SIR'] + TimeSeDT['N_RIS'] + TimeSeDT['N_RIR']
    TimeSeDT['I3'] = TimeSeDT['N_SSI'] + TimeSeDT['N_SRI'] + TimeSeDT['N_RSI'] + TimeSeDT['N_RRI']

    TimeSeDT['R1'] = TimeSeDT['N_RSS'] + TimeSeDT['N_RIS'] + TimeSeDT['N_RSI'] + TimeSeDT['N_RRS'] + TimeSeDT['N_RRI'] + TimeSeDT['N_RSR'] + TimeSeDT['N_RIR'] + TimeSeDT['N_RRR']
    TimeSeDT['R2'] = TimeSeDT['N_SRS'] + TimeSeDT['N_IRS'] + TimeSeDT['N_SRI'] + TimeSeDT['N_RRS'] + TimeSeDT['N_RRI'] + TimeSeDT['N_SRR'] + TimeSeDT['N_IRR'] + TimeSeDT['N_RRR']
    TimeSeDT['R3'] = TimeSeDT['N_SSR'] + TimeSeDT['N_ISR'] + TimeSeDT['N_SIR'] + TimeSeDT['N_RSR'] + TimeSeDT['N_RIR'] + TimeSeDT['N_SRR'] + TimeSeDT['N_IRR'] + TimeSeDT['N_RRR']

    TimeSeDT['Inc1'] = TimeSeDT['Inc_ISS'] + TimeSeDT['Inc_ISR'] + TimeSeDT['Inc_IRS'] + TimeSeDT['Inc_IRR']
    TimeSeDT['Inc2'] = TimeSeDT['Inc_SIS'] + TimeSeDT['Inc_SIR'] + TimeSeDT['Inc_RIS'] + TimeSeDT['Inc_RIR']
    TimeSeDT['Inc3'] = TimeSeDT['Inc_SSI'] + TimeSeDT['Inc_SRI'] + TimeSeDT['Inc_RSI'] + TimeSeDT['Inc_RRI']
    TimeSeDT['IncT'] = TimeSeDT['Inc1'] + TimeSeDT['Inc2'] + TimeSeDT['Inc3']
    TimeSeDT['freq2']=TimeSeDT['Inc2']/TimeSeDT['IncT']

    # EmergeTime = pd.merge(TimeSeDT[TimeSeDT['I2']>0].groupby(['Net','Sigma2','Netrealisation','r2','Sigma','R1t2','deltat','it']).nth(0)[['t','N_SSS','R1']],
    #                                       TimeSeDT[TimeSeDT['I3']>0].groupby(['Net','Sigma2','Netrealisation','r2','Sigma','R1t2','deltat','it']).nth(0)[['t','N_SSS','R1']],
    #                                       on=['Net','Sigma2','Netrealisation','r2','Sigma','R1t2','deltat','it'],suffixes=['_t2','_t3'],how='outer')
    df_em2 = TimeSeDT[TimeSeDT['I2']>0].groupby(['Net','Sigma2','Netrealisation','r2','Sigma','R1t2','deltat','it'], as_index=False).nth(0)
    df_em3 = TimeSeDT[TimeSeDT['I3']>0].groupby(['Net','Sigma2','Netrealisation','r2','Sigma','R1t2','deltat','it'], as_index=False).nth(0)
    EmergeTime = pd.merge(df_em2[['Net','Sigma2','Netrealisation','r2','Sigma','R1t2','deltat','it','t','N_SSS','R1']],
                          df_em3[['Net','Sigma2','Netrealisation','r2','Sigma','R1t2','deltat','it','t','N_SSS','R1']],
                                  on=['Net','Sigma2','Netrealisation','r2','Sigma','R1t2','deltat','it'],suffixes=['_t2','_t3'],how='outer')
    
    EmergeTime.rename({'t_t2':'t2','t_t3':'t3'},axis=1,inplace=True)
    TimeSeDT = pd.merge(TimeSeDT,
                            EmergeTime[['Net','Sigma2','Netrealisation','r2','Sigma','R1t2','deltat','it','t2','t3']], on = ['Net','Sigma2','Netrealisation','r2','Sigma','R1t2','deltat','it'])
    TimeSeDT['t_from_emerge']=TimeSeDT['t']-TimeSeDT['t2']    

    AggonTime_Inc2 = TimeSeDT[TimeSeDT['t_from_emerge']>-2].groupby(['Net','Sigma2','r2','Sigma','R1t2','deltat','t_from_emerge'])[['Inc1','Inc2','freq2']].agg(['count',np.mean,np.std])
    # TimeSeDT[TimeSeDT['t_from_emerge']>-2].groupby(['r2','Sigma','R1t2','deltat','t_from_emerge'])[['Inc1','Inc2','freq2']].describe(percentiles=[0.05,0.50,0.95])
    AggonTime_Inc2.columns = ['_'.join(col) for col in AggonTime_Inc2.columns.values]
    AggonTime_Inc2.reset_index(inplace=True)
    AggonTime_Inc2['inlogit']=np.log(AggonTime_Inc2['freq2_mean']/(1-AggonTime_Inc2['freq2_mean']))
    # AggonTime_Inc2['Inc2_CI']=AggonTime_Inc2['Inc2_95%']-AggonTime_Inc2['Inc2_5%']
    AggonTime_Inc2Max = AggonTime_Inc2.loc[AggonTime_Inc2[AggonTime_Inc2['t_from_emerge']>1].groupby(['Net','Sigma2','r2','Sigma','R1t2','deltat'])['Inc2_mean'].idxmax()]

    TimeIntervalTofitt = pd.merge(AggonTime_Inc2[AggonTime_Inc2['t_from_emerge']==2][['Net','Sigma2','r2','Sigma','R1t2','deltat','Inc2_mean']].rename({'Inc2_mean':'Inc2_init'},axis=1),
             AggonTime_Inc2Max[['Net','Sigma2','r2','Sigma','R1t2','deltat','t_from_emerge','Inc2_mean']].rename({'t_from_emerge':'tmax','Inc2_mean':'Inc2_maxvalue'},axis=1),
            on = ['Net','Sigma2','r2','Sigma','R1t2','deltat'])
    TimeIntervalTofitt['tend']=TimeIntervalTofitt['tmax']
    TimeIntervalTofitt.loc[(TimeIntervalTofitt['tmax']< 1*7), 'tend'] = 1*7

    TimeSeDT = pd.merge(TimeSeDT,TimeIntervalTofitt[['Net','Sigma2','r2','Sigma','R1t2','deltat','tend']], 
                        on = ['Net','Sigma2','r2','Sigma','R1t2','deltat'])

    data = TimeSeDT[(TimeSeDT['IncT']>0)&(TimeSeDT['t_from_emerge']>2)&((TimeSeDT['t_from_emerge']<TimeSeDT['tend']))]
    data.set_index(['Net','Sigma2','r2','Sigma','R1t2','deltat'],inplace=True)
    countzero_Inc = data.groupby(['Net','Sigma2','r2','Sigma','R1t2','deltat'])[['Inc1','Inc2']].agg([('count',len),('countzero' , lambda x: sum(x==0))])
    countzero_Inc.columns = ['_'.join(col) for col in countzero_Inc.columns.values]
    countzero_Inc_indx = countzero_Inc[(countzero_Inc['Inc1_countzero']==countzero_Inc['Inc1_count'])|              (countzero_Inc['Inc2_countzero']==countzero_Inc['Inc2_count'])].index
    data = data.loc[~data.index.isin(countzero_Inc_indx)]
    Models_DT = pd.DataFrame(data.groupby(['Net','Sigma2','r2','Sigma','R1t2','deltat','tend']).apply(model),columns=['model'])
    Models_All = pd.concat([Models_All,Models_DT])

    print('network ',net,' done!', flush=True)
    del TimeSeDT,EmergeTime,data
    gc.collect()

		
Models_All = Models_All.reset_index().set_index(['Net','r2','Sigma2','Sigma','R1t2','deltat','tend'])
coeffResult = Models_All['model'].apply(ModelSumm).stack(level=0).reset_index()
ModelSummDT = Models_All['model'].agg(ModelSumm2).reset_index()

FILE_basename = os.path.basename(FILE)
FILE_out1 = FILE_basename
FILE_out1 = re.sub('Net\d+r\d+',net+'_AllRealisation',FILE_out1)
FILE_out1 = FILE_out1.replace("TimeSerie_","FittedSelCoeff_")
FILE_out1 = re.sub('2Strain.*beta1', '2Strain_beta1', FILE_out1)
FILE_out1 = path+FILE_out1

FILE_out2 = FILE_basename
FILE_out2 = re.sub('Net\d+r\d+',net+'_AllRealisation',FILE_out2)
FILE_out2 = FILE_out2.replace("TimeSerie_","FittedSelCoeffModelSumm_")
FILE_out2 = re.sub('2Strain.*beta1', '2Strain_beta1', FILE_out2)
FILE_out2 = path+FILE_out2

coeffResult.to_csv(FILE_out1, index=False)
ModelSummDT.to_csv(FILE_out2, index=False)

print(FILE_out1, flush=True)



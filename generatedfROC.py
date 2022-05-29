import pandas as pd
from pandas import read_excel

datDir = '/drive3/aliciapm/cnv-project/apm-cnv/tables'
methods = ['ichorcna-cns', 'cnvkit-cns', 'canary-kurtz-arms', 'canary-kurtz-cytobands', 'canary-mse-arms', 'canary-mse-cytobands', 'canary-mse-newcytobands']

dfCNVvalues = read_excel(f'{datDir}/CNV_allcases_r.xlsx')

dfROC = {}
for method in methods:
  Zscore = read_excel(f'{datDir}/{method}-Zscore_table.xlsx')
  merged = pd.merge(dfCNVvalues, Zscore, on='HSTAMP_Label')
  dfROC[method] = merged[["Diagnosis_x", "HSTAMP_Label","t5q_x", "t7q_x", "t11q_x", "t13q_x", "t17p_x", "t20q_x", "t1q_x", "t1p_x", "t14q_x", "ttrisomy8_x", "ttrisomy12_x",
"t5q_y", "t7q_y", "t11q_y", "t13q_y", "t17p_y", "t20q_y", "t1q_y", "t1p_y", "t14q_y", "ttrisomy8_y", "ttrisomy12_y"]]
  dfROC[method].to_excel(f'{datDir}/{method}-dfROC.xlsx')




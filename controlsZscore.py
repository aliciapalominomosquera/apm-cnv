#!/usr/bin/env python3

import pandas as pd
import sys
from CNAdefs import *
from weightedmeanvalue import weightedMeanValues

data = []

# First provide one of the cnv methods to be used: canary-kurtz-cytobands, canary-kurtz-arms, canary-mse-cytobands, canary-mse-newcytobands, canary-mse-arms, wisecondor, testcnvkit, testichor
cnvMethod = 'canary-kurtz-cytobands'
if (len(sys.argv)>=2):
  cnvMethod = sys.argv[1]
startControls = 201
Ncontrols = 50
for i in range(startControls, startControls + Ncontrols + 1):
  HLabel = 'HSTAMP0'+str(i)
  # Default weighted meanvalue (redefine if needed in particular case)
  defaultValue = 'NA'
  # File name and column names used in canary-kurtz output cytobands
  if (cnvMethod=='canary-kurtz-cytobands'):
    fileName = f'/drive3/dkurtz/HEMESTAMP/CANARy/samples/output/Sample_{HLabel}-T1_Tumor.SegmentedGenome.cytobands-noXY.on-off-combined.txt'
    chrColName = "chrNum"; intChromValue = True; startColName = "Start"; endColName = "End"; valueColName = "combinedStoufferZL2CNR"
  # File name and column names used in canary-kurtz output arms
  elif (cnvMethod=='canary-kurtz-arms'):
    fileName = f'/drive3/dkurtz/HEMESTAMP/CANARy/samples/output/Sample_{HLabel}-T1_Tumor.SegmentedGenome.arms.on-off-combined.txt'
    chrColName = "chrNum"; intChromValue = True; startColName = "Start"; endColName = "End"; valueColName = "combinedStoufferZL2CNR" 
    gainThreshold = 1.96; deletionThreshold = -1.96
  # File name and column names used in canary-mse output cytobands
  elif (cnvMethod=='canary-mse-cytobands'):
    fileName = f'canary-python/results-canary/results-canary-mse/Sample_{HLabel}-T1_Tumor.cnvZscores'
    chrColName = "#chr"; intChromValue = False; startColName = "start"; endColName = "end"; valueColName = "gc.corrected.norm.log.std.index.zWeighted.Final"
    gainThreshold = 1.96; deletionThreshold = -1.96
  # File name and column names used in canary-mse new output cytobands
  elif (cnvMethod=='canary-mse-newcytobands'):
    fileName = f'/drive3/mse/CNV/Alicia/results-canary2_new/Sample_{HLabel}-T1_Tumor.cnvZscores'
    chrColName = "#chr"; intChromValue = False; startColName = "start"; endColName = "end"; valueColName = "gc.corrected.norm.log.std.index.zWeighted.Final"
    gainThreshold = 1.96; deletionThreshold = -1.96
  # File name and column names used in canary-mse output arms
  elif (cnvMethod=='canary-mse-arms'):
    fileName = f'/drive3/mse/CNV/Alicia/results-canary5/Sample_{HLabel}-T1_Tumor.cnvZscores'
    chrColName = "#chr"; intChromValue = False; startColName = "start"; endColName = "end"; valueColName = "gc.corrected.norm.log.std.index.zWeighted.Final" 
    gainThreshold = 1.96; deletionThreshold = -1.96
  # File name and column names used in wisecondor output
  elif (cnvMethod=='wisecondor'):
    fileName = f'wisecondor/testSamples/Sample_{HLabel}-T1_Tumor.sorted.samtools-deduped.sorted.offtarget.std.txt'
    chrColName = "chrNum"; intChromValue = True; startColName = "Start"; endColName = "End"; valueColName = "z-score"; defaultValue = 0
    gainThreshold = 1.96; deletionThreshold = -1.96
  # File name and column names used in cnvki-cnst output. Note that cnvkit uses copynumber rather than z-score
  elif (cnvMethod=='cnvkit-cns'):
    fileName = f'cnvkit/results-cnn-tumor/Sample_{HLabel}-T1_Tumor.samtools.call.cns'
    chrColName = "chromosome"; intChromValue = False; startColName = "start"; endColName = "end"; valueColName = "cn"
    gainThreshold = 2.0; deletionThreshold = 2.0
  # File name and column names used in cnvkit-cnr output. Note that cnvkit uses copynumber rather than z-score
  elif (cnvMethod=='cnvkit-cnr'):
    fileName = f'cnvkit/results-cnn-tumor/Sample_{HLabel}-T1_Tumor.samtools.call.cnr'
    chrColName = "chromosome"; intChromValue = False; startColName = "start"; endColName = "end"; valueColName = "cn"
    gainThreshold = 2.0; deletionThreshold = 2.0
 # File name and column names used in ichorcna-cns output. Note that ichorcna uses copynumber rather than z-score
  elif (cnvMethod=='ichorcna-cns'):
    fileName = f'ichorcna/results-ichorcna/{HLabel}.seg'
    chrColName = "chr"; intChromValue = True; startColName = "start"; endColName = "end"; valueColName = "copy.number"
    gainThreshold = 2.0; deletionThreshold = 2.0
 # File name and column names used in ichorcna-cnr output. Note that ichorcna uses copynumber rather than z-score
  elif (cnvMethod=='ichorcna-cnr'):
    fileName = f'ichorcna/results-ichorcna/{HLabel}.cna.seg'
    chrColName = "chr"; intChromValue = True; startColName = "start"; endColName = "end"; valueColName = f'{HLabel}.copy.number'
    gainThreshold = 2.0; deletionThreshold = 2.0
 
  else:
    print(f'Provided cnv method {cnvMethod} not in the supported list: canary-kurtz-cytobands, canary-kurtz-arms, canary-mse-cytobands, canary-mse-newcytobands, canary-mse-arms, wisecondor, cnvkit-cns, cnvkit-cnr, ichorcna-cns,ichorcna-cnr')

  # Obtain zscore for every cna from CNAdefs
  print(f'fileName: {fileName}')
  w_zscores = weightedMeanValues(CNA, cnas, fileName, chrColName, startColName, endColName, valueColName, intChromValue, defaultValue)
  data.append(w_zscores)

dfZscores = pd.DataFrame(data)
#print(dfZscores)
print(dfZscores.mean())
print(dfZscores.std())

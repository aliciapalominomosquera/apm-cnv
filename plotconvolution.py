import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def plotConvolution(value):
  df = pd.read_csv('table-13qconvolution.txt', sep="\t")
  grouped = df.groupby(['width'])
  uid = df.width.unique()
  print(uid)
  #fig, ax = plt.subplots(nrows=len(uid), figsize=(20, 14))
  fig, ax = plt.subplots(figsize=(20, 14))
  for width, df_fit in grouped:
    dfw = df[df['width']==width]
    idmax = dfw[value].idxmax()
    valueMax = dfw[value][idmax]
    centerMax = dfw['center'][idmax]
    intervalMax = dfw['Interval'][idmax]
    print(f"width {width}, max({value}) {valueMax}, interval {intervalMax}, center {centerMax}")
    # get the index of the current width, and use it to index ax
    axi = np.argwhere(uid==width)[0][0]
    # plot to the correct ax based on the index of the width
    ##labels = [ '-'.join(wrap(l, 20)) for l in df_fit.Interval]
    #plot = df_fit.plot(x='center', y=value, ax=ax[axi], label=f'{width}', xlabel='center of the bed file', ylabel=value, title=f'interval width: {width}', marker='.', rot=30) 
    ##plot.set_xticklabels(df_fit.Interval, rotation=20, ha='center', rotation_mode='anchor', fontsize=8)
    ##plot.tick_params(axis="x", pad=20)
    plot = df_fit.plot(x='center', y=value, ax=ax, label=f'w={width}, max {valueMax} in ({intervalMax})', xlabel='center of the bed file', ylabel=value, title=f'{value} distribution - 13q convolution', marker='.', rot=30)
 
    # place the legend outside the plot
    #ax[axi].legend(title='center', bbox_to_anchor=(1.05, 1), loc='upper left')

  plt.tight_layout()
  plt.savefig(f'{value}.png')

plotConvolution('Specificity')
plotConvolution('Sensitivity')
plotConvolution('Accuracy')


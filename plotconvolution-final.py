import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def plotConvolution(value):
  df = pd.read_csv('table-13qconvolution-0.2.txt', sep="\t")
  grouped = df.groupby(['width'])
  uid = df.width.unique()
  print(uid)
  #fig, ax = plt.subplots(nrows=len(uid), figsize=(20, 14))
  fig, ax = plt.subplots(figsize=(10, 7))
  widths = [555423., 1805124.75, 3054826.5, 4304528.25]
  for width in widths:
    dfw = df[df['width']==width]
    idmax = dfw[value].idxmax()
    valueMax = dfw[value][idmax]
    centerMax = dfw['center'][idmax]
    intervalMax = dfw['Interval'][idmax]
    print(f"width {width}, max({value}) {valueMax}, interval {intervalMax}, center {centerMax}")
    plot = dfw.plot(x='center', y=value, ax=ax, label=f'width={width/1000:.0f} kbp', marker='.', rot=30, linewidth=3)  
    plt.xticks(size = 15, rotation = 30 )
    plt.title(f'{value} distribution - 13q convolution', size = 15)
    plt.yticks(size = 15)
    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.xlabel('bed file location', size = 15)
    plt.ylabel(f'{value}', size = 15)
    plt.legend(loc='best', fontsize = 15)
    plt.ylim((0,1))
  plt.tight_layout()
  plt.savefig(f'{value}.png')

plotConvolution('Specificity')
plotConvolution('Sensitivity')
plotConvolution('Accuracy')


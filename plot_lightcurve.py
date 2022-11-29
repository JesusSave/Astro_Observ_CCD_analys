#===================================================================
#
#                      PLOTTING LIGHT CURVES
#
# This script is used to plot light curves from photometry ouput data.
#
# VERSION: 24 Sep 2020
# AUTHOER: QIANG CHEN chen@camk.edu.pl 
#
# PYTHON ENVIRONMENT REQUIREMENT:
# - 
#
# REFERENCE:
# - 
#===================================================================

import os
import numpy as np
import matplotlib.pyplot as plt

PATH = '/work/chuck/chen/obs'

folder,name,deblending = '20190822','stars',False
folder,name,deblending = '20190823','ap',False
folder,name,deblending = '20190827','ap',False
folder,name,deblending = '20190822','ap',False

filters = ['V','U','B','I']

try:
  print('deleting old data')
  os.system(f'rm {PATH}/{folder}/reduced/sorted_{name}_*.png')
except:
  pass

for filt in filters:
  print('SORTING STARS BY EVERY FILTER TYPE:',filt) 
  for i in range(200):
    try:
      f = open(f'{PATH}/{folder}/reduced/sorted_{name}_{filt}{i}.dat','r')
      data = np.array([[float(data) for data in line.split()] for line in f.readlines()])
      f.close()
    except:
      continue

    row,col = data.shape
    if (row<10):
      continue

    print('  filter',filt,'i',i)

    if (name=='stars'):
      x = data[:,1]
      y = data[:,2]
      Julian = data[:,-3]
      #intJul = [int(x) for x in Julian]
      #intJul = int(np.mean(intJul))
      flux = data[:,-5]
    elif(name=='ap' and deblending):
      x = data[:,4]
      y = data[:,5]
      Julian = data[:,-1]
      #intJul = [int(x) for x in Julian]
      #intJul = int(np.mean(intJul))
      flux = data[:,6]
    elif(name=='ap'):
      x = data[:,1]
      y = data[:,2]
      Julian = data[:,-3]
      #intJul = [int(x) for x in Julian]
      #intJul = int(np.mean(intJul))
      flux = data[:,-4]

    #threshold = np.mean(flux)*0.1
    #flux = np.ma.masked_where(flux < threshold, flux)

    plt.clf()
    rows = 1
    cols = 1
    plt.gcf().set_size_inches(8*cols,4*rows)

    plt.subplot(rows,cols,1,title=f'Data: {folder}, Filter: {filt}, ID: {i}, Pixel Coor: ({round(np.mean(x),2)}, {round(np.mean(y),2)})')
    plt.scatter(Julian-int(Julian[0]), flux)
    plt.ylabel('flux')
    plt.xlabel(f'Julian Date -{int(Julian[0])}')
    #plt.ylim(np.mean(flux)-50,np.mean(flux)+50)
    #plt.ylim(0,100)
    #plt.xlim(0.305, 0.40)
    plt.grid(c='0.5',ls=':')
    plt.tick_params(which='both', direction='in', top='true', right='true')
    plt.tight_layout()
    plt.savefig(f'{PATH}/{folder}/reduced/sorted_{name}_{filt}{i}.png',format='png')
    #plt.show()
 

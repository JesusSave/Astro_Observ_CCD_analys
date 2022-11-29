#===================================================================
#
#                       STARS SORTING
#
# This script is used to sorting stars from the table file of each figure.
#
# WARNING:
# 
# VERSION: 24 Sep 2020
# AUTHOER: QIANG CHEN chen@camk.edu.pl 
#
# PYTHON ENVIRONMENT REQUIREMENT:
# - pip install astropy
# - pip install --no-deps ccdproc
# - pip install photutils
# alternatively can use (CAMK):
#   source /home/chen/python-chen/chen3.6/bin/activate
#
# REFERENCE:
# - COMPARISION STARS https://nbviewer.jupyter.org/gist/dokeeffe/416e214f134d39db696c7fdac261964b
#===================================================================


import os
import glob
import numpy as np

PATH = '/work/chuck/chen/obs'

folder,name,deblending = '20190823','ap',False
folder,name,deblending = '20190827','ap',False
folder,name,deblending = '20190822','stars',False
folder,name,deblending = '20190822','ap',False

offset = 100
len_threshold = 10

filters = ['V','U','B','I']
#filters = ['U']

try:
  print('deleting old data')
  os.system(f'rm {PATH}/{folder}/reduced/sorted_{name}_*.dat')
except:
  pass

refs = []
rows = []
for filt in filters:
  print('SORTING STARS BY EVERY FILTER TYPE:',filt) 
  files = glob.glob(f'{PATH}/{folder}/reduced/{name}_light_*{filt}*.dat')

  print('  reading datas')
  row_ref = 0
  datas = []
  data_ref = np.array([])
  for fil in files:
    f = open(fil,'r')
    next(f) # skip header line
    data = np.array([[float(data) for data in line.split()] for line in f.readlines()])
    f.close()
    datas.append(data)
    row,col = data.shape
    print(f'    {fil} contains {row} stars, time {data[0,-1]}')
    if (row >= row_ref):
      row_ref = row
      data_ref = data
      ref = fil
  print('  sorting datas by time ascending')
  datas = sorted(datas, key = lambda x: x[0,-1])
  if (data_ref.size):
    pass
  else:
    print('  NO RESULTS')
    continue
  
  refs.append(ref)
  rows.append(row_ref)
  print(f'  compare to the reference star, search in files within tolerence. Type: {name}')
  ii = 0
  for i in range(row_ref):
    if(name=='ap' and deblending):
      onestar = []
      xcentroid = data_ref[i,4]
      ycentroid = data_ref[i,5]
      flux_ref = data_ref[i,6]
      for data in datas:
        print(f'    ref star {i+1}/{row_ref}, time {data[0,-1]}')
        row,col = data.shape
        for j in range(row):
          x = data[j,4]
          y = data[j,5]
          flux = data[j,6]
          if (xcentroid-offset<x<xcentroid+offset) and (ycentroid-offset<y<ycentroid+offset) and (0.1<flux/flux_ref<10) and (xcentroid<=x):  
            onestar.append(data[j,:])
    elif(name=='ap'):
      onestar = []
      xcentroid = data_ref[i,1]
      ycentroid = data_ref[i,2]
      flux_ref = data_ref[i,3]
      for data in datas:
        print(f'    ref star {i+1}/{row_ref}, time {data[0,-1]}')
        row,col = data.shape
        for j in range(row):
          x = data[j,1]
          y = data[j,2]
          flux = data[j,3]
          if (xcentroid-offset<x<xcentroid+offset) and (ycentroid-offset<y<ycentroid+offset) and (1./10<flux/flux_ref<10):  
            onestar.append(data[j,:])
    elif(name=='stars'):
      onestar = []
      xcentroid = data_ref[i,1]
      ycentroid = data_ref[i,2]
      flux_ref = data_ref[i,-3]
      mag_ref = data_ref[i,-2]
      for data in datas:
        print(f'    ref star {i+1}/{row_ref}, time {data[0,-1]}')
        row,col = data.shape
        for j in range(row):
          x = data[j,1]
          y = data[j,2]
          flux = data[j,-3]
          mag = data[j,-2]
          if (xcentroid-offset<x<xcentroid+offset) and (ycentroid-offset<y<ycentroid+offset) and (1./20<flux/flux_ref<20) and (1./3<mag/mag_ref<3) and (xcentroid<=x):  
            onestar.append(data[j,:])

    if(len(onestar)<len_threshold):
      print('      lack data, abort')
      continue
    else:
      ii+=1
      onestar = np.array(onestar)
      #print('  sort data as time ascending')
      #onestar = onestar[onestar[:,-1].argsort()]
      np.savetxt(f'{PATH}/{folder}/reduced/sorted_{name}_{filt}{ii}.dat', onestar, fmt='%f ', newline='\n')


print('Reference figures with most stars found:')
for k in range(len(refs)):
  ref = refs[k]
  num = rows[k]
  print('  in figure',ref,'found',num,'stars')






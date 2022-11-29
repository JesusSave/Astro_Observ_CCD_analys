#===================================================================
#
#                       STARS SORTING
#
# This script is used to sorting stars from the table file of each figure.
#
# PROCESS:
# 1) Choose Reference Frame
#    select the best quality frame (contains most star) as a reference,
#    pick up several most luminous star (lum_ref_star_limit) as coordinates referencelum_ref_star_limit
# 2) Align Coordinates
#    compare each frame to the reference stars with the nearist triangle pattern,
#    according to the pattern derive transformation matrix, align all the coordinates in that figure
# 3) Sorting Stars
#    iteratively select the star from the reference frame, search in all the other frames within coordinates tolerence
#
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
# - SEARCHING STAR TRIANGLE PATTERNS https://stackoverflow.com/q/43126580/13720066
# - ALIGN POINTS by Homogeneous transformation matrix https://astroalign.readthedocs.io/en/latest/
#===================================================================


def getTriangles(set_X, X_combs):
    """
    Inefficient way of obtaining the lengths of each triangle's side.
    Normalized so that the minimum length is 1.
    """
    triang = []
    for p0, p1, p2 in X_combs:
        d1 = np.sqrt((set_X[p0][0] - set_X[p1][0]) ** 2 +
                     (set_X[p0][1] - set_X[p1][1]) ** 2)
        d2 = np.sqrt((set_X[p0][0] - set_X[p2][0]) ** 2 +
                     (set_X[p0][1] - set_X[p2][1]) ** 2)
        d3 = np.sqrt((set_X[p1][0] - set_X[p2][0]) ** 2 +
                     (set_X[p1][1] - set_X[p2][1]) ** 2)
        d_min = min(d1, d2, d3)
        d_unsort = [d1 / d_min, d2 / d_min, d3 / d_min]
        triang.append(sorted(d_unsort))
    return triang

def sumTriangles(ref_triang, in_triang):
    """
    For each normalized triangle in ref, compare with each normalized triangle
    in B. find the differences between their sides, sum their absolute values,
    and select the two triangles with the smallest sum of absolute differences.
    """
    tr_sum, tr_idx = [], []
    for i, ref_tr in enumerate(ref_triang):
        for j, in_tr in enumerate(in_triang):
            # Absolute value of lengths differences.
            tr_diff = abs(np.array(ref_tr) - np.array(in_tr))
            # Sum the differences
            tr_sum.append(sum(tr_diff))
            tr_idx.append([i, j])
    # Index of the triangles in ref and in with the smallest sum of absolute
    # length differences.
    tr_idx_min = tr_idx[tr_sum.index(min(tr_sum))]
    ref_idx, in_idx = tr_idx_min[0], tr_idx_min[1]
    return ref_idx, in_idx, min(tr_sum)



import os
import glob
import numpy as np
import itertools
import numpy as np
import matplotlib.pyplot as plt
import astroalign as aa

PATH = '/work/chuck/chen/obs'

ALIGN_COOR = False
ALIGN_COOR = True

folder,name,deblending,plot = '20190822','stars',False,False
folder,name,deblending,plot = '20190823','ap',False,False
folder,name,deblending,plot = '20190827','ap',False,False
folder,name,deblending,plot = '20190822','ap',False,False

offset = 10 # better >= 10
len_threshold = 10
lum_ref_star_limit = 10
# >=3, better 5-10, bad figure quality needs bigger, too large in case of slow
in_star_limit = 3 # >=3
flux_ratio = 15 # 
align_accuracy_limit = 5e-3
filters = ['V','B','U','I']


try:
  print('deleting old data')
  os.system(f'rm {PATH}/{folder}/reduced/sorted_{name}_*.dat')
  if(ALIGN_COOR): 
    print('delete align')
    os.system(f'rm {PATH}/{folder}/reduced/aligned_{name}_light_*.dat')
except:
  pass



refs = []
rows = []
datas = {}
align_fails = []
total_files = 0
for filt in filters:
  print('SORTING STARS BY EVERY FILTER TYPE:',filt) 

  print('Reading datas')
  if(ALIGN_COOR):
    files = glob.glob(f'{PATH}/{folder}/reduced/{name}_light_*{filt}*.dat')
  else:
    files = glob.glob(f'{PATH}/{folder}/reduced/aligned_{name}_light_*{filt}*.dat')
  total_files+=len(files)

  for fil in files:
    f = open(fil,'r')
    if(ALIGN_COOR):
      next(f) # skip header line
    data = np.array([[float(data) for data in line.split()] for line in f.readlines()])
    f.close()
    datas[fil] = data

  print('sorting datas by time ascending')
  datas = {k: v for k, v in sorted(datas.items(), key=lambda item: item[1][0,4])}

  print('Setting the most star figure as a reference')
  row_ref = 0
  for key, value in datas.items() :
    row,col = value.shape
    #print(f'    {key} contains {row} stars, time {data[0,4]}')
    if (row > row_ref):
      row_ref = row
      ref = key

  print('sorting reference frame as flux descending')
  if (datas[ref].size):
    if(name=='ap' and deblending):
      datas[ref] = datas[ref][datas[ref][:,6].argsort()[::-1]]
    else:
      datas[ref] = datas[ref][datas[ref][:,3].argsort()[::-1]]
  else:
    print('no reference, reject!')
    continue
  refs.append(ref)
  rows.append(row_ref)



  if(ALIGN_COOR):
    print('Aligning Coordinates')
    for k in range(len(files)):
      k = k
      fil = files[k]
      data = datas[fil]
      row,col = data.shape
      if(row<in_star_limit):
        break
      print(f'  aligning coordinates: {k+1}/{len(files)}\n    reference: {ref}\n    input: {fil}')
      if(name=='ap' and deblending):
        xcentroid = datas[ref][:,4]
        ycentroid = datas[ref][:,5]
        x = data[:,4]
        y = data[:,5]
      else:
        xcentroid = datas[ref][:,1]
        ycentroid = datas[ref][:,2]
        x = data[:,1]
        y = data[:,2]

      set_ref = np.array(list(zip(xcentroid[:lum_ref_star_limit],ycentroid[:lum_ref_star_limit])))
      set_in = np.array(list(zip(x,y)))

      # All possible triangles.
      ref_combs = list(itertools.combinations(range(len(set_ref)), 3))
      in_combs = list(itertools.combinations(range(len(set_in)), 3))

      # Obtain normalized triangles.
      ref_triang, in_triang = getTriangles(set_ref, ref_combs), getTriangles(set_in, in_combs)

      # Index of the ref and in triangles with the smallest difference.
      ref_idx, in_idx, accuracy = sumTriangles(ref_triang, in_triang)
      #print("    smallest difference, accuracy: {}".format(accuracy))
      if (accuracy>align_accuracy_limit):
        print('    quality too low, reject')
        align_fails.append(fil)
        continue

      # Indexes of points in ref and in of the best match triangles.
      ref_idx_pts, in_idx_pts = ref_combs[ref_idx], in_combs[in_idx]
      print ('    triangle ref %s matches triangle in %s' % (ref_idx_pts, in_idx_pts))

      print ("    ref:", [set_ref[_] for _ in ref_idx_pts])
      print ("    input:", [set_in[_] for _ in in_idx_pts])

      ref_pts = np.array([set_ref[_] for _ in ref_idx_pts])
      in_pts = np.array([set_in[_] for _ in in_idx_pts])

      transf, (in_list,ref_list) = aa.find_transform(in_pts, ref_pts)
      transf_in = transf(set_in)
      print(f'    {transf}')    
      if plot:
        plt.scatter(set_ref[:,0],set_ref[:,1], s=100,marker='.', c='r',label='Reference')
        plt.scatter(set_in[:,0],set_in[:,1], s=100,marker='^', c='b',label='Input')
        plt.scatter(transf_in[:,0],transf_in[:,1], s=200,marker='+', c='b',label='Input Aligned')
        plt.plot(ref_pts[:,0],ref_pts[:,1], c='r')
        plt.plot(in_pts[:,0],in_pts[:,1], c='b')
        plt.legend()
        plt.show()
      datas[fil] = np.append(data,transf_in,axis=1)
      head, tail = os.path.split(fil)
      np.savetxt(f'{head}/aligned_{tail}',datas[fil], fmt='%f ', newline='\n')



  print('Sorting Stars')
  print(f'  compare to the reference star, search in files within tolerence. Type: {name}')
  ii = 0
  for i in range(row_ref):
    onestar = []
    print(f'    star {i+1}/{row_ref} in reference: {ref}')
    if(name=='ap' and deblending):
      xcentroid = datas[ref][i,4]
      ycentroid = datas[ref][i,5]
      flux_ref = datas[ref][i,6]
    else:
      xcentroid = datas[ref][i,1]
      ycentroid = datas[ref][i,2]
      flux_ref = datas[ref][i,3]
    for k in range(len(files)):
      fil = files[k]
      data = datas[fil]
      row,col = data.shape
      print(f'    sorting: {k+1}/{len(files)}\n      reference: {ref}\n      ori input: {fil} stars {row}')
      for j in range(row):
        if(name=='ap' and deblending):
          x = data[j,-2]
          y = data[j,-1]
          flux = data[j,6]
        else:
          x = data[j,-2]
          y = data[j,-1]
          flux = data[j,3]
        if (xcentroid-offset<x<xcentroid+offset) and (ycentroid-offset<y<ycentroid+offset): # and (1./flux_ratio<flux/flux_ref<flux_ratio):  
          onestar.append(data[j,:])
    if(len(onestar)<len_threshold):
      print('      lack data, abort')
      continue
    else:
      ii+=1
      onestar = np.array(onestar)
      #print('  sort data as time ascending')
      onestar = onestar[onestar[:,-3].argsort()]
      #print(onestar)
      np.savetxt(f'{PATH}/{folder}/reduced/sorted_{name}_{filt}{ii}.dat', onestar, fmt='%f ', newline='\n')



print('Reference figures with most stars found:')
for k in range(len(refs)):
  ref = refs[k]
  num = rows[k]
  print('  in figure',ref,'found',num,'stars')

if(ALIGN_COOR):
  for fail in align_fails:
    print('  align fails in', fail)
  print(f'align fail total number {len(align_fails)} / {total_files}')

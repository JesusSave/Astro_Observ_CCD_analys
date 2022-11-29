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
folder = '20190823'
folder = '20190827'
folder = '20190822'


star22 = [['stars', 'B', '10', 'b.'],
        ['stars', 'V', '15', 'g.'],
       ] #22 

star27 = [['stars', 'B', '36', 'b.'],
        ['stars', 'V', '47', 'g.'],
       ] #27

ap22 =  [['ap', 'B', '6', 'b.'],
       ['ap', 'V', '6', 'g.'],
      ] #22

ap23 =  [['ap', 'B', '3', 'b.'],
       ['ap', 'V', '4', 'g.'],
      ] #23

ap27 =  [['ap', 'B', '3', 'b.'],
       ['ap', 'V', '4', 'g.'],
      ] #27

cases = star27
cases = ap27

plt.clf()
rows = 1
cols = 1
plt.gcf().set_size_inches(8*cols,4*rows)
plt.subplot(rows,cols,1)

mx = []
my = []
for k in range(len(cases)):
    name, filt, ID, fmt =  cases[k]

    f = open(f'{PATH}/{folder}/reduced/sorted_{name}_{filt}{ID}.dat','r')
    data = np.array([[float(data) for data in line.split()] for line in f.readlines()])
    f.close()
    x = data[:,1]
    y = data[:,2]
    mx.append(np.mean(x))
    my.append(np.mean(y))
    Julian = data[:,-3]
    intJul = [int(x) for x in Julian]
    intJul = int(np.mean(intJul))
    if (name=='stars'):
      flux = data[:,-3]
    else:
      flux = data[:,-4]
    #threshold = np.mean(flux)*0.1
    #flux = np.ma.masked_where(flux < threshold, flux)
    #if (np.mean(flux)<0):
    #  flux = flux + max(abs(flux))

    plt.errorbar(Julian-intJul, flux/max(flux),fmt = fmt,label=f'Filter: {filt}')
    #plt.errorbar(Julian-intJul, flux/max(flux), yerr = np.std(flux)/max(flux) ,fmt = fmt,label=f'Filter: {filt}')
    #plt.errorbar(Julian-intJul, flux/max(flux), yerr = np.sqrt(flux)/max(flux) ,fmt = fmt,label=f'Filter: {filt}')

plt.title(f'Observation date: {folder}, AB variable at pixels ({round(np.mean(mx),2)}, {round(np.mean(my),2)})')
plt.ylabel('normalized flux')
plt.xlabel(f'Julian Date -{intJul}')
#plt.ylim(np.mean(flux)-50,np.mean(flux)+50)
#plt.ylim(0,100)
#plt.xlim(0.305, 0.40)
plt.grid(c='0.5',ls=':')
plt.legend(loc=4)
plt.tick_params(which='both', direction='in', top='true', right='true')
plt.tight_layout()
plt.savefig(f'{PATH}/{folder}/reduced/final_AB_{name}_{folder}.png',format='png')
plt.show()
 

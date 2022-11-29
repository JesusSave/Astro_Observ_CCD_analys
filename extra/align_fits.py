'''
import drizzlepac
from drizzlepac import tweakreg

#unlearn tweakreg
#unlearn imagefindpars
tweakreg.TweakReg('*.fit',threshold=100,searchrad=1e6, tolerance=10)
#tweakreg.TweakReg('*.fit',threshold=10)#,searchrad=4.0,\
'''

import drizzlepac
from drizzlepac import tweakreg
tweakreg.TweakReg('*.fit',
       tolerance=15, # matches
       searchrad=10, # 1 rad ~ 20 pixel 
       #residplot='None', # do not plot residual
       interactive=False, # do not show or ask !
       imagefindcfg={'threshold' : 200, 'conv_width' : 3.5},
       refimagefindcfg={'threshold' : 400, 'conv_width' : 2.5},
       updatehdr=True, shiftfile=True, outshifts='shift.txt')

'''
import drizzlepac
from drizzlepac import tweakreg
tweakreg.TweakReg('*.fit',
       imagefindcfg=dict(threshold=10, conv_width=3.5),
       refimagefindcfg=dict(threshold=10, conv_width=2.5),
       updatehdr=False, shiftfile=True, outshifts='shift.txt')
'''





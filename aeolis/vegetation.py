'''This file is part of AeoLiS.
   
AeoLiS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
   
AeoLiS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with AeoLiS.  If not, see <http://www.gnu.org/licenses/>.
   
AeoLiS  Copyright (C) 2015 Bas Hoonhout

bas.hoonhout@deltares.nl         b.m.hoonhout@tudelft.nl
Deltares                         Delft University of Technology
Unit of Hydraulic Engineering    Faculty of Civil Engineering and Geosciences
Boussinesqweg 1                  Stevinweg 1
2629 HVDelft                     2628CN Delft
The Netherlands                  The Netherlands

'''
from __future__ import absolute_import, division

import logging
import numpy as np

# package modules
import aeolis.wind
#from aeolis.utils import *

# initialize logger
logger = logging.getLogger(__name__)

def initialize (s,p):

    s['rhoveg'][:, :] = p['veg_file']

    if np.isnan(s['rhoveg'][0, 0]):
        s['rhoveg'][:,:] = 0.

    ix = s['rhoveg'] < 0.03
    s['rhoveg'][ix] *= 0.
    s['hveg'][:, :] = p['hveg_max']*np.sqrt(s['rhoveg'])

    s['germinate'][:,:] = (s['rhoveg']>0)
    s['lateral'][:,:] = 0.

    return s

def germinate (s,p):
    
    # time
    
    n = (365.25*24*3600/p['dt'])
    
    # Germination
    
    p_germinate_year = p['germinate'] # 0.05 / year 
    p_germinate_dt = 1-(1-p_germinate_year)**(1/n)
    germination = np.random.random((s['germinate'].shape))
    
    s['germinate'] += (s['dzveg'] >= -0.01) * (p['_time'] > p['dz_interval']) * (germination <= p_germinate_dt)
    s['germinate'] = np.minimum(s['germinate'],1.)

    # SINGLE EMBRYO TEST
#    s['germinate'][25,25] = 1.

    # Lateral expension
    
    dx = s['dsu'][2,2]
    
    p_lateral_year = p['lateral']  #0.2  / year
    p_lateral_dt = 1-(1-p_lateral_year)**(1/n)
    p_lateral_cell = 1 - (1-p_lateral_dt)**(1/dx)
    
    drhoveg = np.zeros((p['ny']+1, p['nx']+1, 4))
    
    drhoveg[:,1:,0] = np.maximum((s['rhoveg'][:,:-1]-s['rhoveg'][:,1:]) / s['dsu'][:,1:], 0.)  # positive x-direction
    drhoveg[:,:-1,1] = np.maximum((s['rhoveg'][:,1:]-s['rhoveg'][:,:-1]) / s['dsu'][:,:-1], 0.)  # negative x-direction
    drhoveg[1:,:,2] = np.maximum((s['rhoveg'][:-1,:]-s['rhoveg'][1:,:]) / s['dnu'][1:,:], 0.)  # positive y-direction
    drhoveg[:-1,:,3] = np.maximum((s['rhoveg'][1:,:]-s['rhoveg'][:-1,:]) / s['dnu'][:-1,:], 0.)  # negative y-direction
    
    lat_veg = drhoveg > 0.
    
    s['dxrhoveg'] = np.sum(lat_veg[:,:,:],2)
    
    p_lateral = p_lateral_cell * s['dxrhoveg']
    
    s['lateral'] += (germination <= p_lateral)
    s['lateral'] = np.minimum(s['lateral'],1.)

    return s

def grow (s, p):
    
    ix = np.logical_or(s['germinate'] != 0., s['lateral'] != 0.) * ( p['V_ver'] > 0.)

    s['drhoveg'][:,:] *= 0.

    # Reduction of vegetation growth due to sediment burial
    s['dhveg'][ix] = p['V_ver'] * (1 - s['hveg'][ix]/p['hveg_max']) - np.abs(s['dzveg'][ix])*p['veg_gamma'] # m/year
    s['dhveg'] = np.maximum(s['dhveg'], -3.)

    # if p['_time'] > p['dz_interval']:
    #     import matplotlib.pyplot as plt
    #     plt.pcolormesh(s['x'], s['y'], s['dhveg'], vmin=-1, vmax=1)
    #     plt.colorbar()
    #     plt.show()

    # Adding growth
    s['hveg'] += s['dhveg']*p['dt']/(365.25*24*3600)

    # Compute the density

    s['hveg'] = np.maximum(np.minimum(s['hveg'], p['hveg_max']), 0.)
    s['rhoveg'] = (s['hveg']/p['hveg_max'])**2

    # jx = np.logical_and(s['lateral'] != 0.,s['rhoveg'] < 0.02)
    # s['rhoveg'][jx] += 0.4*s['dxrhoveg'][jx]*0.5

    # Plot has to vegetate again after dying

    s['germinate'] *= (s['rhoveg']!=0.)
    s['lateral'] *= (s['rhoveg']!=0.)

    # Dying of vegetation due to hydrodynamics (Dynamic Vegetation Limit)

    s['rhoveg'] *= (s['zb']+0.01 >= s['zs'])
    s['hveg'] *= (s['zb'] + 0.01 >= s['zs'])
    s['germinate'] *= (s['zb']+0.01 >=s['zs'])
    s['lateral'] *= (s['zb']+0.01 >=s['zs'])
    
    return s

def vegshear(s, p):

    # Raupach, 1993
    
    roughness = 16.
    nf = p['nfractions']

# -----------------------------
    ets = np.zeros(s['zb'].shape)
    etn = np.zeros(s['zb'].shape)
    ets[:,:] = 1.

    ix = s['ustar'] != 0.
    ets[ix] = s['ustars'][ix]/s['ustar'][ix]
    etn[ix] = s['ustarn'][ix]/s['ustar'][ix]

    s['vegfac'] = 1. / np.sqrt(1. + roughness * s['rhoveg'])

    s['ustar'] *= s['vegfac']
    s['ustars']  = s['ustar'] * ets
    s['ustarn']  = s['ustar'] * etn
# ------------------------------
    
    # s['vegfac']= 1./(1. + roughness*s['rhoveg'])
    #
    # s['tau'] *= s['vegfac']
    #
    # s['taus'] = ets * s['tau']
    # s['taun'] = etn * s['tau']
    #
    # s['ustar'] = np.sqrt(s['tau'] /p['rhoa'])
    # s['ustars']  = s['ustar'] * ets
    # s['ustarn']  = s['ustar'] * etn

# --------------------------------------
    
    # vegfac = np.repeat(s['vegfac'][:,:,np.newaxis], nf, axis = 0)
    # tauth = p['rhoa']*s['uth']**2
    # tauth /= vegfac
    #
    # s['uth'] = np.sqrt(tauth/p['rhoa'])
    
    
    return s
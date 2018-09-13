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
    
    s['hveg'][:,:] = 0.
    s['rhoveg'][:,:] = 0.
    s['germinate'][:,:] = 0.
    s['lateral'][:,:] = 0.
    
    return s

def germinate (s,p):
    
    # time
    
    n = (30*24*3600/p['dt'])
#    n_year = (365.25*24*3600/p['dt'])
    
    # Germination
    
    p_germinate_month = p['germinate'] # 0.005 / month 
    p_germinate_dt = 1-(1-p_germinate_month)**(1/n)
    germination = np.random.random((s['germinate'].shape))
    
    s['germinate'] += (s['dz_avg'] >= 0.) * (p['_time'] > p['dz_interval']) * (germination <= p_germinate_dt)
    
    s['germinate'][:10,:] = 0.
    s['germinate'][-10:,:] = 0.
    s['germinate'][:,:10] = 0.
    s['germinate'][:,-10:] = 0.
    
    s['germinate'] = np.minimum(s['germinate'],1.)
    
    # Lateral expension
    
    p_lateral_month = 0.5 #0.99  / month
    p_lateral_dt = 1-(1-p_lateral_month)**(1/n)
    
    drhoveg = np.zeros((p['ny']+1, p['nx']+1, 4))
    
    drhoveg[:,1:,0] = np.maximum((s['rhoveg'][:,:-1]-s['rhoveg'][:,1:]) / s['dsu'][:,1:], 0.)  # positive x-direction
    drhoveg[:,:-1,1] = np.maximum((s['rhoveg'][:,1:]-s['rhoveg'][:,:-1]) / s['dsu'][:,:-1], 0.)  # negative x-direction
    drhoveg[1:,:,2] = np.maximum((s['rhoveg'][:-1,:]-s['rhoveg'][1:,:]) / s['dnu'][1:,:], 0.)  # positive y-direction
    drhoveg[:-1,:,3] = np.maximum((s['rhoveg'][1:,:]-s['rhoveg'][:-1,:]) / s['dnu'][:-1,:], 0.)  # negative y-direction

    s['dxrhoveg'] = (np.sum(drhoveg[:,:,:],2)>0.)
    
    #Determine slopes
    
    slope = np.zeros((p['ny']+1, p['nx']+1, 4))
    
    slope[:,1:,0] = np.abs((s['zb'][:,:-1]-s['zb'][:,1:]) / s['dsu'][:,1:])  # positive x-direction
    slope[:,:-1,1] = np.abs((s['zb'][:,1:]-s['zb'][:,:-1]) / s['dsu'][:,:-1])  # negative x-direction
    slope[1:,:,2] = np.abs((s['zb'][:-1,:]-s['zb'][1:,:]) / s['dnu'][1:,:])  # positive y-direction
    slope[:-1,:,3] = np.abs((s['zb'][1:,:]-s['zb'][:-1,:]) / s['dnu'][:-1,:])  # negative y-direction
    
    limit_slope = 15. # 15
    limit_slope_ratio = np.arctan(np.deg2rad(limit_slope))
    
    slope_total = np.sum(slope[:,:,:],2)
    
    p_lateral = p_lateral_dt * s['dxrhoveg'] * (slope_total <= limit_slope_ratio)
    s['lateral'] += (germination <= p_lateral)
    s['lateral'] = np.minimum(s['lateral'],1.)
    
    s['lateral'][:10,:] = 0.
    s['lateral'][-10:,:] = 0.
    s['lateral'][:,:10] = 0.
    s['lateral'][:,-10:] = 0.

    return s

def grow (s, p):
    
    ix = np.logical_or(s['germinate'] != 0., s['lateral'] != 0.)
    
    V_ver = p['V_ver']/(365.25*24*3600)
#    V_lat = p['V_lat']/(365.25*24*3600)

    s['drhoveg'][:,:] *= 0.
    
    # Determine vegetation cover gradients
    
#    drhoveg = np.zeros((p['ny']+1, p['nx']+1, 4))
#    
#    drhoveg[:,1:,0] = np.maximum((s['rhoveg'][:,:-1]-s['rhoveg'][:,1:]) / s['dsu'][:,1:], 0.)  # positive x-direction
#    drhoveg[:,:-1,1] = np.maximum((s['rhoveg'][:,1:]-s['rhoveg'][:,:-1]) / s['dsu'][:,:-1], 0.)  # negative x-direction
#    drhoveg[1:,:,2] = np.maximum((s['rhoveg'][:-1,:]-s['rhoveg'][1:,:]) / s['dnu'][1:,:], 0.)  # positive y-direction
#    drhoveg[:-1,:,3] = np.maximum((s['rhoveg'][1:,:]-s['rhoveg'][:-1,:]) / s['dnu'][:-1,:], 0.)  # negative y-direction

    # Reduction of vegetation growth due to sediment burial
    short_factor = 0.8
    dz_opt  = p['dz_opt']   #[m/year]
    dz_tol  = p['dz_tol']   #[m/year]
    dz      = s['dzyear_avg']*short_factor   #[m/year]
    
    s['sedfac']=np.maximum((1-((dz-dz_opt)/dz_tol)**2),-2.)
    
    # Growth of vegetation
    s['drhoveg'][ix]    = V_ver*s['sedfac'][ix]
    
    # Growth towards limit (upper and lower)
    
    ix = s['drhoveg'] > 0.
    s['drhoveg'][ix] *= (1-s['rhoveg'][ix]) 
    ix = s['drhoveg'] < 0.
    s['drhoveg'][ix] *= s['rhoveg'][ix]
    
#    s['drhoveg']        += V_lat*np.sum(drhoveg[:,:,:],2)
    
    s['rhoveg'] += s['drhoveg']*p['dt']
    s['rhoveg'] = np.maximum(np.minimum(s['rhoveg'],1.),0.)
    
    s['germinate'] *= (s['rhoveg']!=0.)
    s['lateral'] *= (s['rhoveg']!=0.)
    
    # Dying of vegetation
    
    s['rhoveg'] *= (s['zb']+0.1 >= s['zs'])
    s['germinate'] *= (s['zb']+0.1 >=s ['zs'])
    s['lateral'] *= (s['zb']+0.1 >=s ['zs'])
    
    return s

def vegshear(s, p):
    
    roughness = 16.
    
    s['vegfac']= 1./(1. + roughness*s['rhoveg'])
    
    #filter
#    s = aeolis.wind.filter_low(s, p, 'vegfac', 'x', 2.)
#    s = aeolis.wind.filter_low(s, p, 'vegfac', 'y', 2.)
#    
    ets = np.zeros(s['zb'].shape)
    etn = np.zeros(s['zb'].shape)
    ets[:,:] = 1.
#    
    ix = s['ustar'] != 0.
    ets[ix] = s['ustars'][ix]/s['ustar'][ix]
    etn[ix] = s['ustarn'][ix]/s['ustar'][ix]
    
    s['ustar'] *= s['vegfac']
#    s['ustarn'] *= s['vegfac']
    
    s['ustars'] = ets * s['ustar']
    s['ustarn'] = etn * s['ustar']
    
#    s['ustar'] = np.hypot(s['ustars'], s['ustarn'])
    
    return s

#def grow_cdm (s, p, l, t):
# 
# CDM Method 1
#    Vveg_year = 15. # m/year
#    Vveg = Vveg_year / (365.25 * 24 * 3600) # m/s
#    Hveg = 1. # [m]
#    
#    # From Marco 2011
#    
#    Vveg *= (s['dz_avg'] > 0.)
#    s['dhveg'] = Vveg * (1- s['hveg']/Hveg) - np.abs(s['dz_avg'])/p['dt']
#    s['hveg'] += s['dhveg']*p['dt']
#    
#    s['hveg'] = np.maximum(s['hveg'],0.)
#
#    s['rhoveg'] = (s['hveg']/Hveg)**2.
#    
#    s['drhoveg'] = (1-s['rhoveg'])/(p['tveg']*3600*24) - p['eroveg']/p['Hveg']*s['dz_avg'] # Make month or something...
#    s['rhoveg'] += s['drhoveg']*p['dt']
#    
#    # Germination has to happen again after vegetation died
#    s['germinate'] *= (s['rhoveg']!=0.)
#    
#    return s
#    
#def grow_cdm2 (s, p):
#    
#    ix = s['germinate'] != 0.
#    
#    s['drhoveg'][:,:] *= 0.
#    s['drhoveg'][ix] = (1.-s['rhoveg'][ix])/(p['tveg']*86400.) - p['eroveg']*np.abs(s['dz_avg'][ix])/(p['Hveg']*p['dt']) # Make month or something... 
#    
#    s['rhoveg'] += s['drhoveg']*p['dt']
#    s['rhoveg'] = np.maximum(np.minimum(s['rhoveg'],1.),0.)
#    
#    s['germinate'] *= (s['rhoveg']!=0.)
#    
#    return s
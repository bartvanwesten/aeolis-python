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
#from aeolis.utils import *

# initialize logger
logger = logging.getLogger(__name__)

def germinate (s,p,l,t):
    
    s['dzdt'] = s['zb'] - s['zbold']
    s['germinate'] += (t>0)*(s['zb']>=s['zbold']) #18000000
    s['germinate'] = np.minimum(s['germinate'],1.)
    
    return s
        
def grow (s, p, l, t):
    
    tveg = 50.*24.*3600. # [s]
    gamma = 1. # [-]
    Hveg = 1. # [m]

    s['dvegrho'] = (((1.-s['vegrho'])/tveg)-(gamma/Hveg)*np.abs((s['zb']-s['zbold']))/p['dt'])*s['germinate']
    
    s['vegrho'] += s['dvegrho']*p['dt']
    s['vegrho'] = np.minimum(s['vegrho'],1.)
    s['vegrho'] = np.maximum(s['vegrho'],0.)
    
#    s['germinate'] *= (s['vegrho']==0.)
    
    return s
    
def vegshear(s, p):
    
    roughness = 16. #16.
    
    s['vegfac']= 1./(1. + roughness*s['vegrho'])
    s['tau'] *= s['vegfac']
    
#    m_kCut_veg = 1. #20
#    
#    taufft = np.fft.fft2(s['tau'])
#    dk = 2.0 * np.pi / (np.max(s['x']))
#    taufft *= np.exp(-(dk*s['x'])**2./(2.*m_kCut_veg**2.))
#    s['tau'] = np.maximum(np.real(np.fft.ifft2(taufft)),0)
    
    return s
    
def initialize (s,p):
    
    s['vegrho'][:,:] = 0.
    s['germinate'][:,:] = 0.
    
    return s
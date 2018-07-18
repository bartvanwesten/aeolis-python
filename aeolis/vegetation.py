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
    
    s['germinate'] += (s['dz_avg'] > 0.)
    s['germinate'] = np.minimum(s['germinate'],1.)
    
    return s
        
def grow_cdm (s, p, l, t):
    
    Vveg_year = 15. # m/year
    Vveg = Vveg_year / (365.25 * 24 * 3600) # m/s
    Hveg = 1. # [m]
    
    # From Marco 2011
    
    Vveg *= (s['dz_avg'] > 0.)
    s['dhveg'] = Vveg * (1- s['hveg']/Hveg) - np.abs(s['dz_avg'])/p['dt']
    s['hveg'] += s['dhveg']*p['dt']
    
    s['hveg'] = np.maximum(s['hveg'],0.)

    s['rhoveg'] = (s['hveg']/Hveg)**2.
    
    # Germination has to happen again after vegetation died
    s['germinate'] *= (s['rhoveg']!=0.)
    
    return s
    
def grow (s, p):
    
    
    
    return s

def vegshear(s, p):
    
    roughness = 16.
    
    s['vegfac']= 1./(1. + roughness*s['rhoveg'])
    
    ets = s['ustars']/s['ustar']
    etn = s['ustarn']/s['ustar']
    
    s['ustar'] *= s['vegfac']
    
    s['ustars'] = ets * s['ustar']
    s['ustarn'] = etn * s['ustar']
    
    return s
    
def initialize (s,p):
    
    s['hveg'][:,:] = 0.
    s['rhoveg'][:,:] = 0.
    s['germinate'][:,:] = 0.
    
    return s
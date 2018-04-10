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
from aeolis.utils import *


# initialize logger
logger = logging.getLogger(__name__)

def initialize(s, p):

    ny = p['ny']

    # initialize x-dimensions
    s['x'][:,:] = p['xgrid_file']
    
    # World coordinates of z-points
    s['sz'][:,:] = s['x'][:,:] # Different from XBeach
    
    # World coordinates of u-points
    s['su'][:,1:-1] = 0.5 * (s['sz'][:,:-1] + s['sz'][:,1:])
    s['su'][:,0]    = 1.5 *  s['sz'][:,0]   - 0.5 * s['sz'][:,1]
    s['su'][:,-1]   = 1.5 *  s['sz'][:,-1]  - 0.5 * s['sz'][:,-2]
    
    # World coordinates of v-points
    s['sv'][1:-1,:] = 0.5 * (s['sz'][:-1,:] + s['sz'][1:,:])
    s['sv'][0,:]    = 1.5 *  s['sz'][0,:]   - 0.5 * s['sz'][1,:]
    s['sv'][-1,:]   = 1.5 *  s['sz'][-1,:]  - 0.5 * s['sz'][-2,:]
    
    # World coordinates of c-points
    s['sc'][1:-1,1:-1] = 0.25 *(s['sz'][:-1,:-1] + s['sz'][:-1,1:] + s['sz'][1:,:-1] + s['sz'][1:,1:])
    s['sc'][1:-1,0]    = 0.5 * (s['su'][:-1,0]  + s['su'][1:,0])
    s['sc'][1:-1,-1]   = 0.5 * (s['su'][:-1,-1] + s['su'][1:,-1])
    s['sc'][0,1:-1]    = 0.5 * (s['sv'][0,:-1]  + s['sv'][0,1:])
    s['sc'][-1,1:-1]   = 0.5 * (s['sv'][-1,:-1] + s['sv'][-1,1:])
    
    s['sc'][0,0]   = s['su'][0,0]   # Different from XBeach
    s['sc'][0,-1]  = s['su'][0,-1]  # Different from XBeach
    s['sc'][-1,0]  = s['su'][-1,0]  # Different from XBeach
    s['sc'][-1,-1] = s['su'][-1,-1] # Different from XBeach
    
    # Distances
    s['dsz'][:,:] = ((s['su'][:,:-1]-s['su'][:,1:])**2.+(s['nu'][:,:-1]-s['nu'][:,1:])**2.)**0.5
    s['dsu'][:,:] = ((s['sz'][:,:-1]-s['sz'][:,1:])**2.+(s['nz'][:,:-1]-s['nz'][:,1:])**2.)**0.5
    
    s['dsv'][:,:] = ((s['sc'][:,:-1]-s['sc'][:,1:])**2.+(s['nc'][:,:-1]-s['nc'][:,1:])**2.)**0.5
    s['dsc'][:,:] = ((s['sv'][:,:-1]-s['sv'][:,1:])**2.+(s['nv'][:,:-1]-s['nv'][:,1:])**2.)**0.5
    
    # initialize y-dimension
    if ny == 0:
        s['y'][:,:] = 0.
        s['nz'][:,:] = 0.
        s['nu'][:,:] = 0.
        s['nv'][:,:] = 0.
        s['dnz'][:,:] = 1.
        s['dnu'][:,:] = 1.
        s['dnv'][:,:] = 1.
        s['dnc'][:,:] = 1.
        s['alfaz'][:,:] = 0.
    else:
        # initialize y-dimensions
        s['y'][:,:] = p['ygrid_file']
        
        # World coordinates of z-points
        s['nz'][:,:] = s['y'][:,:] # Different from XBeach
        
        # World coordinates of u-points
        s['nu'][:,1:-1] = 0.5 * (s['nz'][:,:-1] + s['nz'][:,1:])
        s['nu'][:,0]    = 1.5 *  s['nz'][:,0]   - 0.5 * s['nz'][:,1]
        s['nu'][:,-1]   = 1.5 *  s['nz'][:,-1]  - 0.5 * s['nz'][:,-2]
        
        # World coordinates of v-points
        s['nv'][1:-1,:] = 0.5 * (s['nz'][:-1,:] + s['nz'][1:,:])
        s['nv'][0,:]    = 1.5 *  s['nz'][0,:]   - 0.5 * s['nz'][1,:]
        s['nv'][-1,:]   = 1.5 *  s['nz'][-1,:]  - 0.5 * s['nz'][-2,:]
        
        # World coordinates of c-points
        s['nc'][1:-1,1:-1] = 0.25 *(s['nz'][:-1,:-1] + s['nz'][:-1,1:] + s['nz'][1:,:-1] + s['nz'][1:,1:])
        s['nc'][0,1:-1]    = 0.5 * (s['nv'][0,:-1]  + s['nv'][0,1:])
        s['nc'][-1,1:-1]   = 0.5 * (s['nv'][-1,:-1] + s['nv'][-1,1:])
        s['nc'][1:-1,0]    = 0.5 * (s['nu'][:-1,0]  + s['nu'][1:,0])
        s['nc'][1:-1,-1]   = 0.5 * (s['nu'][:-1,-1] + s['nu'][1:,-1])
        
        s['nc'][0,0]   = s['nv'][0,0]   # Different from XBeach
        s['nc'][0,-1]  = s['nv'][0,-1]  # Different from XBeach
        s['nc'][-1,0]  = s['nv'][-1,0]  # Different from XBeach
        s['nc'][-1,-1] = s['nv'][-1,-1] # Different from XBeach
        
        # Distances
        s['dnz'][:,:] = ((s['nv'][:-1,:]-s['nv'][1:,:])**2.+(s['sv'][:-1,:]-s['sv'][1:,:])**2.)**0.5
        s['dnu'][:,:] = ((s['sc'][:-1,:]-s['sc'][1:,:])**2.+(s['nc'][:-1,:]-s['nc'][1:,:])**2.)**0.5
        s['dnv'][:,:] = ((s['sz'][:-1,:]-s['sz'][1:,:])**2.+(s['nz'][:-1,:]-s['nz'][1:,:])**2.)**0.5
        s['dnc'][:,:] = ((s['su'][:-1,:]-s['su'][1:,:])**2.+(s['nu'][:-1,:]-s['nu'][1:,:])**2.)**0.5

    # Cell areas
    s['dsdnu'][:,:] = (0.5*(s['dsc'][:-1,:]+s['dsc'][1:,:])) * (0.5*(s['dnz'][:,:-1]+s['dnz'][:,1:]))
    s['dsdnv'][:,:] = (0.5*(s['dsz'][:-1,:]+s['dsz'][1:,:])) * (0.5*(s['dnc'][:,:-1]+s['dnc'][:,1:]))
    s['dsdnz'][:,:] = (0.5*(s['dsv'][:-1,:]+s['dsv'][1:,:])) * (0.5*(s['dnu'][:,:-1]+s['dnu'][:,1:]))
    
    # Inverse cell areas
    s['dsdnui'][:,:] = 1. / s['dsdnu']
    s['dsdnvi'][:,:] = 1. / s['dsdnv']
    s['dsdnzi'][:,:] = 1. / s['dsdnz']
    
    # Alfaz, grid orientation in z-points
    s['alfaz'][1:-1,:] = np.arctan2(s['x'][2:,:] - s['x'][:-2,:], s['y'][2:,:] - s['y'][:-2,:])
    s['alfaz'][0,:] = s['alfaz'][1,:]
    s['alfaz'][-1,:] = s['alfaz'][-2,:]
    
    print(s['sz'][:,:])
    print(s['nz'][:,:])
#    print(s['sv'][:,:])
#    print(s['sc'][:,:])
#    print(s['dsz'][:,:])
#    print(s['dsu'][:,:])
#    print(s['dsv'][:,:])
#    print(s['dsc'][:,:])
    print(s['dsdnz'][:,:])
    print(s['dsdnu'][:,:])
    
    return s
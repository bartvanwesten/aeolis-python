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

'''
The separation bubble is modified from the CDM source code.

Missing parts:
    - Smoothening of bubble
    - Frequency filter in y-direction
    - Cutting of separation bubble of new dune, only one separation bubble for each slice in y-direction.
'''
def separation_shear(s,p):
      
    #Calculate delta for relation separation bubble and shear stress
    
    nx = p['nx'] + 1
    ny = p['ny'] + 1
    x = s['x']
    dx = x[0,1]-x[0,0]
    nf = p['nfractions']
    
    # ACCORDING CDM
    m_tau_sepbub = .05
    slope = np.tan(np.deg2rad(p['Mcr_dyn']))*dx
    delta = 1./(slope*m_tau_sepbub)
    
    dzsep = np.maximum(s['zsep'] - s['zb'], 0.)
    
    s['zsepdelta'] = np.minimum(np.maximum(1. - delta * dzsep, 0.),1.)
    
    if p['process_separation']:
        s['ustar'] *=s['zsepdelta']
        zsepdelta = np.repeat(s['zsepdelta'][:,:,np.newaxis], nf, axis = 0)
        s['uu'] *= zsepdelta   
    
        # Calculte ustar and tau
        
        ets = np.zeros(s['zb'].shape)
        etn = np.zeros(s['zb'].shape)
            
        ix = s['ustar'] != 0.
        ets[ix] = s['ustars'][ix]/s['ustar'][ix]
        etn[ix] = s['ustarn'][ix]/s['ustar'][ix]
            
        s['tau'] = p['rhoa']*s['ustar']**2
    
        s['taus'] = s['tau']*ets
        s['taun'] = s['tau']*etn
        
        s['ustars'] = s['ustar']*ets
        s['ustarn'] = s['ustar']*etn
        
        # uu
        
        ets = np.zeros(s['uu'].shape)
        etn = np.zeros(s['uu'].shape)
            
        ix = s['uu'] != 0.
        ets[ix] = s['uus'][ix]/s['uu'][ix]
        etn[ix] = s['uun'][ix]/s['uu'][ix]
        
        s['uus'] = s['uu']*ets
        s['uun'] = s['uu']*etn
        
#    s['Ts'] *= s['zsepdelta']
#    s['Ts'] = np.maximum(s['Ts'],0.001)
        
#    s['tau'] = np.hypot(s['taus'],s['taun'])
#    s['ustar'] = np.hypot(s['ustars'],s['ustarn'])
    
    return s

def separation(s, p):
    
    if p['process_separation'] :
        
        # Initialize grid and bed dimensions
         
        nx = p['nx'] + 1
        ny = p['ny'] + 1
        x = s['x']
        z = s['zb']
        dx = x[0,1]-x[0,0]
        
        # Initialize arrays
    
        dz  = np.zeros((ny, nx))
        stall = np.zeros((ny, nx))
        bubble = np.zeros((ny, nx)) 
        zfft = np.zeros((ny, nx), dtype=np.complex)
        k = np.array(range(0,nx))
        
        # Compute angle of separation
        
        dz[:,:-1] = np.rad2deg(np.arctan((z[:,1:]-z[:,:-1])/dx))
        dz[:,-1] = dz[:,-2]
        
        # Determine location of separation bubbles
        stall += np.logical_and(dz < 0, abs(dz) > p['M_sep'])
        stall[1:-1,:] += np.logical_and(stall[1:-1,:]==0, stall[:-2,:]>0, stall[2:,:]>0)
        stall[:,1:-1] += np.logical_and(stall[:,1:-1]==0, stall[:,:-2]>0, stall[:,2:]>0)  
        
        for j in range(0,ny): 
            for i in range(0,nx):
                if (np.sum(stall[j,i-1])==0. and stall[j,i]==1):
                    bubble[j, i-1] = 1. # -0
                    
#         In order to reduce the amount of separation bubbles in y-direction to one. JUST TEMP!
#            for i in range(0,nx):
#                if np.sum(bubble[j,:i-1])>0.:
#                    bubble[j, i] = 0.
                
        s['stall']=stall
        s['bubble']=bubble
        
        s['zsep'][:,:] = np.min(s['zb'])
        s['zsepnofil'][:,:] = np.min(s['zb'])
        
        zmin = np.zeros(ny)
        
        # Walk through all separation bubbles and determine polynoms
        
        for j in range(0,ny):
            
            zmin[j] = np.min(z[j,:])
            
            for i in range(0,nx):
                    
                if bubble[j,i]==1:
                    
                    h = s['zb'] - zmin[j]
                    
                    xb = x[j,i]
                    hb = h[j,i]
                    
                    dhdx0 = (z[j,i]-z[j,i-1])/(1*dx)
                    
                    # Separation bubble for shear stresses
                    
                    s['zsepshear'][j,:] = zmin[j]
                    h, k_max = poly(s, p, i, j, hb, xb, dx, nx, dhdx0, x, h, 11.)
                    s['zsepshear'][j,:] += (h[j,:]+zmin[j])
                    
                    # Zero order polynom
                    
                    h, k_max = poly(s, p, i, j, hb, xb, dx, nx, dhdx0, x, h, p['M_dSlope'])
                    s['zsepnofil'][j,:] += h[j,:]
                    
                    # Filter high frequencies (NEW FILTER REQUIRED)
                            
                    dk = 2.0 * np.pi / (np.max(s['x']))
                    zfft[j,:] = np.fft.fft(h[j,:])
                    zfft[j,:] *= np.exp(-(dk*k*dx)**2./(2.*p['m_kCut']**2.))
                    h[j,:] = np.real(np.fft.ifft(zfft[j,:]))
                    
                    # First order polynom
                    
                    dhdx1 = (h[j,i]-h[j,i-1])/dx
                    h, k_max = poly(s, p, i, j, hb, xb, dx, nx, dhdx1, x, h, p['M_dSlope'])
                    s['zsep'][j,i:k_max] = np.maximum(h[j,i:k_max],s['zsep'][j,i:k_max])
                    s['zsep'][j,:] += zmin[j]
        
        s['dzsep'] = np.maximum(s['zsep'] - s['zb'],0.)

    return s

def poly(s, p, i, j, hb, xb, dx, nx, dhdx, x, h, slope):
    
    slope = np.deg2rad(slope)
    
    a = dhdx / slope
    
    l = np.minimum(np.maximum((1.5 * hb / slope) * (1 + a*0.25 + 0.125*a**2),.1),200.)
    
    a3 = 2*hb/l**3 + dhdx/l**2
    a2 = -3*hb/l**2 - 2*dhdx/l
    
    k_max = min(i+int(l/dx),int(nx))
    xs = x[j,i:k_max] - xb
    
    h[j,i:k_max] = ((a3*xs + a2) * xs + dhdx)*xs + hb
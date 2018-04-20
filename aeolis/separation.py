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

def separation(s, p):
    
    # Initialize grid and bed dimensions
     
    nx = p['nx'] + 1
    ny = p['ny'] + 1
    x = s['x']
    z = s['zb']
    dx = x[0,1]-x[0,0]
    zmin = np.min(z)
    h = s['zb'] - zmin
    
    # Calculate delta for relation separation bubble and shear stress
    
    m_tau_sepbub = 0.05
    slope = np.tan(np.deg2rad(p['Mcr_dyn']))*dx
    delta = 1./(slope*m_tau_sepbub)
    
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
    
    bubble[:,:-3] += np.logical_and(stall[:,2:-1]==1, stall[:,1:-2]==0, z[:,2:-1]>z[:,3:])
    bubble[:,1:-2] = bubble[:,3:]
    bubble[:,-2:] = 0
    
    # Walk through all separation bubbles and determine polynoms
    
    for j in range(0,ny):
        for i in range(0,nx):
                
            if bubble[j,i]==1:
                
                xb = x[j,i]
                hb = h[j,i]
                
                # Zero order polynom

                dhdx = (z[j,i]-z[j,i-2])/(2.*dx)
                    
                h = poly(s, p, i, j, hb, xb, dx, nx, dhdx, x, h)
                
                # Filter high frequencies
                        
                dk = 2.0 * np.pi / (np.max(s['x']))
                zfft[j,:] = np.fft.fft(h[j,:])
                zfft[j,:] *= np.exp(-(dk*k*dx)**2./(2.*p['m_kCut']**2.))
                h[j,:] = np.real(np.fft.ifft(zfft[j,:]))
                
                # First order polynom
                
                s['zsep'][j,:] = 0.
                
                dhdx = (h[j,i]-h[j,i-1])/dx

                h = poly(s, p, i, j, hb, xb, dx, nx, dhdx, x, h)
                
                #
    
    s['zsep'] += zmin
    s['zsepdelta'] = np.minimum(np.maximum(1. - delta * (s['zsep'] - s['zb']),
                     np.zeros((ny, nx))),
                     np.zeros((ny, nx)) + 1)
    
    

    return s

def poly(s, p, i, j, hb, xb, dx, nx, dhdx, x, h):
    
    a = dhdx / p['M_dSlope']
    
    l = np.maximum((1.5 * hb / p['M_dSlope']) * (1 + a*0.25 + 0.125*a**2),0.1)
    
    a3 = 2*hb/l**3 + dhdx/l**2
    a2 = -3*hb/l**2 - 2*dhdx/l
    
    k_max = min(i+int(l/dx),int(nx))
    
    xs = x[j,i:k_max] - xb
    
    s['zsep'][j,i:k_max] = ((a3*xs + a2) * xs + dhdx)*xs + hb
    
    h[j,i:k_max] = s['zsep'][j,i:k_max]
        
    return h
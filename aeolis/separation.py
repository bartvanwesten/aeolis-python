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
    m_tau_sepbub = .5
    slope = np.tan(np.deg2rad(p['Mcr_dyn']))*dx
    delta = 1./(slope*m_tau_sepbub)
    
    s['zsepdelta'] = np.minimum(np.maximum(1. - delta * (s['dzsep']), 0.),1.)
    
    if p['process_separation']:
#        s['tau'] *=s['zsepdelta']
        s['taus'] *=s['zsepdelta']
        s['taun'] *=s['zsepdelta']
        s['ustar'] *=s['zsepdelta']
        s['ustars'] *=s['zsepdelta']
        s['ustarn'] *=s['zsepdelta']
        
        zsepdelta = np.repeat(s['zsepdelta'][:,:,np.newaxis], nf, axis = 0)
        
        s['uu'] *= zsepdelta
        s['uus'] *= zsepdelta
        s['uun'] *= zsepdelta
        
#        s['Ts'] *= s['zsepdelta']
#        s['Ts'] = np.maximum(s['Ts'],0.01)
        
    s['tau'] = np.hypot(s['taus'],s['taun'])
    s['ustar'] = np.hypot(s['ustars'],s['ustarn'])
    
    return s

def separation(s, p):
    
    if p['process_separation']:
        
        # Initialize grid and bed dimensions
         
        nx = p['nx'] + 1
        ny = p['ny'] + 1
        x = s['x']
        z = s['zb']
        dx = x[0,1]-x[0,0]
        
        # Initialize arrays
    
        dz  = np.zeros((ny, nx))
#        ddz = np.zeros((ny, nx))
        stall = np.zeros((ny, nx))
        bubble = np.zeros((ny, nx)) 
    #    zsepcrit = np.zeros((ny, nx))
        zfft = np.zeros((ny, nx), dtype=np.complex)
        k = np.array(range(0,nx))
        
        # Compute angle of separation
        
        dz[:,:-1] = np.rad2deg(np.arctan((z[:,1:]-z[:,:-1])/dx))
        dz[:,-1] = dz[:,-2]
        
#        ddz[:,:-1] = dz[:,1:]-dz[:,:-1]
#        ddz[:,-1] = ddz[:,-2] 
        
        s['dz']=dz
        
        # Determine location of separation bubbles
    #    if p['separationcrit'] == 'down':
        stall += np.logical_and(dz < 0, abs(dz) > p['M_sep'])
        stall[1:-1,:] += np.logical_and(stall[1:-1,:]==0, stall[:-2,:]>0, stall[2:,:]>0)
        stall[:,1:-1] += np.logical_and(stall[:,1:-1]==0, stall[:,:-2]>0, stall[:,2:]>0)  
        
        for j in range(0,ny): 
            for i in range(0,nx):
                if (np.sum(stall[j,i-1])==0. and stall[j,i]==1):
                    bubble[j, i-2] = 1.
                    
#         In order to reduce the amount of separation bubbles in y-direction to one. JUST TEMP!
            for i in range(0,nx):
                if np.sum(bubble[j,:i-1])>0.:
                    bubble[j, i] = 0.
                
        s['stall']=stall
        s['bubble']=bubble
        
        s['zsep'][:,:] = 0.
        s['zsepnofil'][:,:] = 0.
        
        zmin = np.zeros(ny)
        
        # Walk through all separation bubbles and determine polynoms
        
        for j in range(0,ny):
            
            zmin[j] = np.min(z[j,:])
            
            for i in range(0,nx):
                    
                if bubble[j,i]==1:
                    
                    h = s['zb'] - zmin[j]
                    
                    xb = x[j,i]
                    hb = h[j,i]
                    
                    # Zero order polynom
    
                    dhdx0 = (z[j,i]-z[j,i-2])/(2.*dx)
                    
#                    s['zsepnofil'][j,:] = 0.
                        
                    h, k_max = poly(s, p, i, j, hb, xb, dx, nx, dhdx0, x, h, p['M_dSlope'])
                    
#                    s['zsepnofil'][j,:] += h[j,:]
                    
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
                    
                    # Area limit
#                    area = np.sum(s['zsep'][j,i:k_max]) - np.sum(s['zb'][j,i:k_max])
#                    if area < 2.:
#                        s['zsep'][j,:] = zmin[j]
                    
                    # Separation bubble for shear stresses
                    
                    s['zsepshear'][j,:] = 0.
                    h = poly2(s, p, i, j, hb, xb, dx, nx, dhdx0, x, h, 0.8*p['M_dSlope'])
                    s['zsepshear'][j,:] += zmin[j]                    
        
        s['dzsep'] = np.maximum(s['zsep'] - s['zb'],0.)

    return s

def poly(s, p, i, j, hb, xb, dx, nx, dhdx, x, h, slope):
    
    slope = np.deg2rad(slope)
    
    a = dhdx / slope
    
    l = np.maximum((1.5 * hb / slope) * (1 + a*0.25 + 0.125*a**2),.1)
    
    a3 = 2*hb/l**3 + dhdx/l**2
    a2 = -3*hb/l**2 - 2*dhdx/l
    
    k_max = min(i+int(l/dx),int(nx))
    
    xs = x[j,i:k_max] - xb
    
#    if l==2.0:
#        s['zsep'][j,i:k_max]=0.
#    else:
    h[j,i:k_max] = ((a3*xs + a2) * xs + dhdx)*xs + hb
    
#    h[j,i:k_max] = s['zsep'][j,i:k_max]
        
    return h, k_max

def poly2(s, p, i, j, hb, xb, dx, nx, dhdx, x, h, slope):
    
    a = dhdx / slope
    
    l = np.maximum((1.5 * hb / slope) * (1 + a*0.25 + 0.125*a**2),.1)
    a3 = 2*hb/l**3 + dhdx/l**2
    a2 = -3*hb/l**2 - 2*dhdx/l
    
    k_max = min(i+int(l/dx),int(nx))
    
    xs = x[j,i:k_max] - xb
    
    s['zsepshear'][j,i:k_max] = ((a3*xs + a2) * xs + dhdx)*xs + hb
    
    h[j,i:k_max] = s['zsepshear'][j,i:k_max]
        
    return h
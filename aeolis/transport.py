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
    
def grainspeed(s,p):
    '''Compute velocity threshold and grain speed according to Duran 2007 (p. 42)

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters

    Returns
    -------
    dict
        Spatial grids
        '''
  
    # Create each grain fraction
    
    nf = p['nfractions']
    d = p['grain_size']
    A = 0.95
    B = 5.12
    srho = p['rhop']/p['rhoa']
    
    g       = np.repeat(p['g'], nf, axis = 0)
    v       = np.repeat(p['v'], nf, axis = 0)
    alfa    = np.repeat(p['alfa'], nf, axis = 0)
    A       = np.repeat(A, nf, axis = 0)
    B       = np.repeat(B, nf, axis = 0)
    kappa   = np.repeat(p['karman'], nf, axis = 0)
    rhop    = np.repeat(p['rhop'], nf, axis = 0)
    rhoa    = np.repeat(p['rhoa'], nf, axis = 0)
    
    srho    = rhop / rhoa
    
    uth     = np.zeros(d.shape)
    uth     = p['A'] * np.sqrt((p['rhop'] - p['rhoa']) / p['rhoa'] * g * d)
    
    # ustar (no vector) is defined in Parteli 2013 as the average shear velocity (True?) (Not spatially varying?)
    # For now it is assumed to be the shear velocity over a flat bed
    ustar0  = np.average(s['ustar0'])
    ustar0  = np.repeat(ustar0[np.newaxis], nf, axis = 0)
    
    ustar   = np.repeat(s['ustar'][:,:,np.newaxis], nf, axis = 0)
    ustars  = np.repeat(s['ustars'][:,:,np.newaxis], nf, axis = 0)
    ustarn  = np.repeat(s['ustarn'][:,:,np.newaxis], nf, axis = 0)

    # Drag coefficient (Duran, 2007 -> Jimenez and Madsen, 2003)
    
    Sstar   = np.zeros(d.shape)
    Cd      = np.zeros(d.shape)
    uf      = np.zeros(d.shape)
    
    Sstar   = d/(4*v)*np.sqrt(g*d*(srho-1.))
    Cd      = (4/3)*(A+np.sqrt(2*alfa)*B/Sstar) 
    
    uf = np.sqrt(4/(3*Cd)*(srho-1)*g*d)
    
    # Compute saltation heights
    
    r       = 1. # Duran 2007, p. 33
    c       = 14./(1.+1.4*r)
    c       = np.repeat(c, nf, axis = 0)
    
    tv      = (v/g**2)**(1/3) # +- 5.38

    zm      = c * uth * tv  # characteristic height of the saltation layer +- 20 mm
    z0      = d/20. # grain based roughness layer +- 10 mu m
    
    z1      = 0.003 # reference height
    z1      = np.repeat(z1, nf, axis = 0)
    
    # Efficient wind velocity

    ueff = np.zeros(s['uth'].shape)

#    ueff = (uth/kappa)*(np.log(z1/z0)+2*(np.sqrt(1+z1/zm*(ustar**2/uth**2-1))-1))
    ueff = (uth/kappa)*(np.log(z1/z0)+ z1/zm*(ustar/uth-1))
    ueff0 = (uth/kappa)*(np.log(z1/z0)+ z1/zm*(ustar0/uth-1))
    
    # Grain velocity
    
    dhs = np.zeros(s['uth'].shape)
    dhn = np.zeros(s['uth'].shape)
    ets = np.zeros(s['uth'].shape)
    etn = np.zeros(s['uth'].shape)
    dh  = np.zeros(s['uth'].shape)
    et  = np.zeros(s['uth'].shape)
    A   = np.zeros(s['uth'].shape)
    
    dhs0 = np.zeros(s['uth'].shape)
    dhn0 = np.zeros(s['uth'].shape)
    ets0 = np.zeros(s['uth'].shape)
    etn0 = np.zeros(s['uth'].shape)
    dh0  = np.zeros(s['uth'].shape)
    et0  = np.zeros(s['uth'].shape)
    A0   = np.zeros(s['uth'].shape)
    
    x       = np.repeat(s['x'][:,:,np.newaxis], nf, axis = 0)
    y       = np.repeat(s['y'][:,:,np.newaxis], nf, axis = 0)
    z       = np.repeat(s['zb'][:,:,np.newaxis], nf, axis = 0)
    
    # Improve according to new gridparams! IS this correct?
    
    dhs[:,:-1,:] = (z[:,1:,:]-z[:,:-1,:])/(x[:,1:,:]-x[:,:-1,:])
    dhn[1:-1,:,:] = (z[2:,:,:]-z[:-2,:,:])/(y[2:,:,:]-y[:-2,:,:])
    
    dhs[:,:,:] = 0.
    dhn[:,:,:] = 0.
    
    # Boundaries
    dhs[:,0,:] = dhs[:,1,:]
    dhn[0,:,:] = dhn[1,:,:]    
#    dhs[:,-1,:] = dhs[:,-2,:]
    dhn[-1,:,:] = dhn[-2,:,:]
    
    s['dhs']=dhs
    s['dhn']=dhn
    
    ets[:,:,:] += 1.
    etn[:,:,:] += 0.
    
    ix = ustar != 0.
    
    ets[ix] = ustars[ix] / ustar[ix]
    etn[ix] = ustarn[ix] / ustar[ix]
    
    ets0 = 1.
    etn0 = 0.

    dh[:,:,:] = np.hypot(dhs,dhn)
    et[:,:,:] = np.hypot(ets,etn)
    
    dh0[:,:,:] = 1.
    et0[:,:,:] = 1.
    
    A = et + 2*alfa*dh
    A0 = et0 + 2*alfa*dh0
    
#    As = ets + 2*alfa*dhs
#    An = etn + 2*alfa*dhn
#    
#    A = np.hypot(As,An)

    s['uus'] = (ueff-uf/(np.sqrt(2*alfa)*A))*ets-(np.sqrt(2*alfa)*uf/A)*dhs
    s['uun'] = (ueff-uf/(np.sqrt(2*alfa)*A))*etn-(np.sqrt(2*alfa)*uf/A)*dhn
    
    s['uus0'] = (ueff0-uf/(np.sqrt(2*alfa)*A0))*ets0-(np.sqrt(2*alfa)*uf/A0)*dhs0
    s['uun0'] = (ueff0-uf/(np.sqrt(2*alfa)*A0))*etn0-(np.sqrt(2*alfa)*uf/A0)*dhn0
    
    s['uu'] = np.hypot(s['uus'], s['uun'])
    s['uu0'] = np.hypot(s['uus0'], s['uun0'])

    return s

def equilibrium(s, p):
    '''Compute equilibrium sediment concentration following Bagnold (1937)

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters

    Returns
    -------
    dict
        Spatial grids

    '''

    if p['process_transport']:
        
        nf = p['nfractions']
        ustar = s['ustar'][:,:,np.newaxis].repeat(nf, axis=2)
        ustar0 = s['ustar0'][:,:,np.newaxis].repeat(nf, axis=2)

        ix = ustar != 0.
        
        # equilibrium concentration
        
        s['Cu'] = np.zeros(ustar.shape)
        
        if p['method_transport'].lower() == 'bagnold':
            s['Cu'][ix] = np.maximum(0., p['Cb'] * p['rhoa'] / p['g'] \
                                     * (ustar[ix] - s['uth'][ix])**3 / s['uu'][ix])
            s['Cu0'] = np.maximum(0., p['Cb'] * p['rhoa'] / p['g'] \
                                     * (ustar0 - s['uth0'])**3 / s['uu0'])
        elif p['method_transport'].lower() == 'kawamura':
            s['Cu'][ix] = np.maximum(0., p['Cb'] * p['rhoa'] / p['g'] \
                                     * (ustar[ix] + s['uth'][ix])**2 * (ustar[ix] - s['uth'][ix]) / s['uu'][ix])
        elif p['method_transport'].lower() == 'lettau':
            s['Cu'][ix] = np.maximum(0., p['Cb'] * p['rhoa'] / p['g'] \
                                     * (ustar[ix] - s['uth'][ix]) * ustar[ix]**2 / s['uu'][ix])
        else:
            logger.log_and_raise('Unknown transport formulation [%s]' % p['method_transport'], exc=ValueError)
    
    s['Cu'] *= p['accfac']

    return s


def compute_weights(s, p):
    '''Compute weights for sediment fractions

    Multi-fraction sediment transport needs to weigh the transport of
    each sediment fraction to prevent the sediment transport to
    increase with an increasing number of sediment fractions. The
    weighing is not uniform over all sediment fractions, but depends
    on the sediment availibility in the air and the bed and the bed
    interaction parameter ``bi``.

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters

    Returns
    -------
    numpy.ndarray
        Array with weights for each sediment fraction

    '''

    w_air = normalize(s['Ct'], s['Cu'])
    w_bed = normalize(s['mass'][:,:,0,:], axis=2)

    w = (1. - p['bi']) * w_air \
        + (1. - np.minimum(1., (1. - p['bi']) * np.sum(w_air, axis=2, keepdims=True))) * w_bed
    w = normalize(w, axis=2)
    
    return w, w_air, w_bed


def renormalize_weights(w, ix):
    '''Renormalizes weights for sediment fractions

    Renormalizes weights for sediment fractions such that the sum of
    all weights is unity. To ensure that the erosion of specific
    fractions does not exceed the sediment availibility in the bed,
    the normalization only modifies the weights with index equal or
    larger than ``ix``.

    Parameters
    ----------
    w : numpy.ndarray
        Array with weights for each sediment fraction
    ix : int
        Minimum index to be modified

    Returns
    -------
    numpy.ndarray
        Array with weights for each sediment fraction

    '''
    
    f = np.sum(w[:,:,:ix], axis=2, keepdims=True)
    w[:,:,ix:] = normalize(w[:,:,ix:], axis=2) * (1. - f)

    # normalize in case of supply-limitation
    # use uniform distribution in case of no supply
    w = normalize(w, axis=2, fill=1./w.shape[2])

    return w

def saturation_time(s,p):
    
    # THIS PART HAS TO BE CHECKED (BART) - frac
    
    nf      = p['nfractions']
    d       = p['grain_size']
    rgamma  = .2
    
    g       = np.repeat(p['g'], nf, axis = 0)
    alfa    = np.repeat(p['alfa'], nf, axis = 0)
    rgamma  = np.repeat(rgamma, nf, axis = 0)
    
    uth = s['uth0']
    ustar  = np.repeat(s['ustar'][:,:,np.newaxis], nf, axis = 0)
    
    zb  = np.repeat(s['zb'][:,:,np.newaxis], nf, axis = 0)
    zsep  = np.repeat(s['zsep'][:,:,np.newaxis], nf, axis = 0)
    
    uth0 = np.average(s['uth0'])
    ustar0 = np.average(s['ustar0'])
    
    # Saturation time (Direct implementation from Duran 2007)

    Ts      = np.zeros(ustar.shape)
    Ts      = 2*s['uu']*alfa/(g*rgamma)/((ustar/uth)**2-1)
    
    Tsmax   = 3.
    Tsmin   = 0.01
    
    ix      = ustar <= uth0
    Ts[ix]  = 3.0
    
#    T0      = 1.75 * ustar0 / g   
#    Ts0     = T0 * (uth0/(ustar0-uth0)) / rgamma[0]

    ix      = zsep > 0.1
    Ts[ix]  = .01

    # Range
    Ts = np.maximum(np.minimum(Ts,Tsmax),Tsmin) 
    s['Ts'][:,:] = p['T']#Ts[:,:,0]

    # Filter
#    
#    Cut = 5.0
#    
#    parfft = np.fft.fft2(s['Ts'])
#    dk = 2.0 * np.pi / np.max(s['x'])
#    fac = np.exp(-(dk*s['x']**2.)/(2.*Cut**2.))
#    parfft *= fac
#    s['Ts'] = np.real(np.fft.ifft2(parfft))
#    
#    for j in range(0,p['ny']):
#        tauavg = np.average(s['Ts'][j,:10])
#        s['Ts'][j,:]+=(Ts0-tauavg)
#    for i in range(0,p['nx']):
#        tauavg = np.average(s['Ts'][:10,i])
#        s['Ts'][:,i]+=(Ts0-tauavg)
##        
#    Ts = p['T'] * (uth/(ustar-uth)) / rgamma
#    s['Ts'][:,:] = p['T'] #np.minimum(np.maximum(Ts[:,:,0],0.01),3.)
    
#    print('step')
    
    return s
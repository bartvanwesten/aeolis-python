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
import aeolis.wind

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
#    alfa    = np.repeat(p['alfa'], nf, axis = 0)
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
    
    r       = 1. # Duran 2007, p. 33
    c       = 14./(1.+1.4*r)
    c       = np.repeat(c, nf, axis = 0)
    
    tv      = (v/g**2)**(1/3) # +- 5.38
    
    lv      = (v**2/(p['A']**2*g*(srho-1)))**(1/3)

    zm      = c * uth * tv  # characteristic height of the saltation layer +- 20 mm
    z0      = d/20. # grain based roughness layer +- 10 mu m
    z1      = 35. * lv # reference height
    
    alfa    = 0.17 * d / lv
    p['alfa'] = alfa
    
    Sstar   = d/(4*v)*np.sqrt(g*d*(srho-1.))
    Cd      = (4/3)*(A+np.sqrt(2*alfa)*B/Sstar)**2
    
    uf = np.sqrt(4/(3*Cd)*(srho-1)*g*d)
    
    # Efficient wind velocity

    ueff = np.zeros(s['uth'].shape)

    ueff = (uth/kappa)*(np.log(z1/z0)+ z1/zm*(np.maximum(ustar/uth,1.)-1))
    
    ueff0 = (uth/kappa)*(np.log(z1/z0)+ z1/zm*(ustar0/uth-1))
    
    # Grain velocity
    
    ets = np.zeros(s['uth'].shape)
    etn = np.zeros(s['uth'].shape)
    
    dhs0 = np.zeros(s['uth'].shape)
    dhn0 = np.zeros(s['uth'].shape)
    ets0 = np.zeros(s['uth'].shape)
    etn0 = np.zeros(s['uth'].shape)
    dh0  = np.zeros(s['uth'].shape)
    et0  = np.zeros(s['uth'].shape)
    Ax0   = np.zeros(s['uth'].shape)
    
    dzs = np.zeros(s['zb'].shape)
    dzn = np.zeros(s['zb'].shape)
    
    #slopes
    
    z = s['zb'].copy()
    x = s['x']
    y = s['y']
    
#    ix = s['zsep'] > z
#    z[ix] = s['zsep'][ix]
    
    dzs[:,1:-1] = 0.5*(z[:,2:]-z[:,:-2])/(x[:,2:]-x[:,:-2])
    dzn[1:-1,:] = 0.5*(z[:-2,:]-z[2:,:])/(y[:-2,:]-y[2:,:])
    
    # Boundaries
    dzs[:,0] = dzs[:,1]
    dzn[0,:] = dzn[1,:]    
    dzs[:,-1] = dzs[:,-2]
    dzn[-1,:] = dzn[-2,:]
    
#    dhs   = np.repeat(dzs[:,:,np.newaxis], nf, axis = 0)/(2*np.sqrt(2))
#    dhn   = np.repeat(dzn[:,:,np.newaxis], nf, axis = 0)/(2*np.sqrt(2))
    
    dhs   = np.repeat(dzs[:,:,np.newaxis], nf, axis = 0)
    dhn   = np.repeat(dzn[:,:,np.newaxis], nf, axis = 0)
    
    # Wind direction
    
    ix = ustar != 0.
    
    ets[ix] = ustars[ix] / ustar[ix]
    etn[ix] = ustarn[ix] / ustar[ix]
    
    ets0 = 1.
    etn0 = 0.
    
    # Add
    
    dh0[:,:,:] = 0.
    et0[:,:,:] = 1.
    
    Axs = ets + 2*alfa*dhs
    Axn = etn + 2*alfa*dhn
    Ax = np.hypot(Axs, Axn)
    
    Ax0 = et0 + 2*alfa*dh0
    
    # Compute grain speed
    
    ix = Ax != 0.

    s['uus'][ix] = (ueff[ix]-uf/(np.sqrt(2*alfa)*Ax[ix]))*ets[ix]-(np.sqrt(2*alfa)*uf/Ax[ix])*dhs[ix]
    s['uun'][ix] = (ueff[ix]-uf/(np.sqrt(2*alfa)*Ax[ix]))*etn[ix]-(np.sqrt(2*alfa)*uf/Ax[ix])*dhn[ix]
    
    s['uus'] = np.maximum(s['uus'] , 0.)
    
    s['uus0'] = (ueff0-uf/(np.sqrt(2*alfa)*Ax0))*ets0-(np.sqrt(2*alfa)*uf/Ax0)*dhs0
    s['uun0'] = (ueff0-uf/(np.sqrt(2*alfa)*Ax0))*etn0-(np.sqrt(2*alfa)*uf/Ax0)*dhn0
    
    s['uu'] = np.hypot(s['uus'], s['uun'])
    s['uu0'] = np.hypot(s['uus0'], s['uun0'])
    
##    ix = s['zsepdelta'] < 0.9
#    s['uus'][ix,0] *= 0.
#    s['uun'][ix,0] *= 0.
#    s['uus'][ix,0] *= 0.
    
    s['uus'][:,:,0] *= s['zsepdelta']
    s['uun'][:,:,0] *= s['zsepdelta']
    s['uu'][:,:,0] *= s['zsepdelta']

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
        elif p['method_transport'].lower() == 'cdm':
            s['Cu'][ix] = np.maximum(0., p['alfa'] * p['rhoa'] / p['g'] \
                                     * (ustar[ix]**2 - s['uth'][ix]**2))
            s['Cu0'][ix] = np.maximum(0., p['alfa'] * p['rhoa'] / p['g'] \
                                     * (ustar0[ix]**2 - s['uth0'][ix]**2))
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
    rgamma  = .2
    
    tauth = p['rhoa']*s['uth0']**2
    tau  = np.repeat(s['tau'][:,:,np.newaxis], nf, axis = 0)
    tau0  = np.repeat(s['tau0'][:,:,np.newaxis], nf, axis = 0)
    ustar  = np.repeat(s['ustar'][:,:,np.newaxis], nf, axis = 0)
#    ustar0  = np.repeat(s['ustar0'][:,:,np.newaxis], nf, axis = 0)
    
    zsepdelta = np.repeat(s['zsepdelta'][:,:,np.newaxis], nf, axis = 0)
    
    # Saturation time (Direct implementation from Duran 2007)

    Ts = 2*s['uu']*p['alfa']/(p['g'])*np.maximum((tauth/((tau-tauth)*rgamma)),0.)
    s['Ts'] = Ts[:,:,0]
    
    Ts_max = 2.
    Ts_min = 0.05
    
    uth = p['A'] * np.sqrt((p['rhop'] - p['rhoa']) / p['rhoa'] * p['g'] * p['grain_size'][0])
    
    # Low grain speeds
    ix  = s['ustar'] <= uth
    s['Ts'][ix]  = Ts_max
    
    jx = s['zsep'] > s['zb']
    s['Ts'][jx] = 0.25
    
    s['Ts'] = np.maximum(np.minimum(s['Ts'],Ts_max),Ts_min)
    
    return s
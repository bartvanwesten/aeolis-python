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

    nf = p['nfractions']
    
    #From shear to shear velocity
    
    ustar = s['ustar']
    ustars = s['ustars']
    ustarn = s['ustarn']
    
#    uth0 = np.zeros(s['uth0'].shape)
    uth0 = s['uth0'][:,:,0]
    
    # Collect variables # uf and d has to be improved! (BART)
    
    d = p['grain_size'][0]
    srho = p['rhop']/p['rhoa']
    g = p['g']
    v = p['v']
    alfa = p['alfa']
    A = 0.95
    B = 5.12
    kappa = p['karman']

    # Drag coefficient (Duran, 2007 -> Jimenez and Madsen, 2003)
    
#    Sstar = np.zeros(d.shape)
#    Cd = np.zeros(d.shape)
    
    Sstar = d/(4*v)*np.sqrt(g*d*(srho-1.))
    Cd = (4/3)*(A+np.sqrt(2*alfa)*B/Sstar) 
    
    # Grain settling velocity

#    uf = np.zeros(d.shape)
    uf = np.sqrt(4/(3*Cd)*(srho-1)*g*d)
    
    # Efficient wind velocity

    r = 1. # Duran 2007, p. 33
    c = 14./(1.+1.4*r)
    tv = (v/g**2)**(1/3) # +- 5.38
    
    zm = np.zeros(ustar.shape)
#    z0 = np.zeros(d.shape)
    
    zm = c * uth0 * tv  # characteristic height of the saltation layer +- 20 mm
    z0 = d/20. # grain based roughness layer +- 10 mu m
    z1 = 0.003 # reference height
    
    ueff = np.zeros(ustar.shape)
    ueff = (uth0/kappa)*(np.log(z1/z0)+z1/zm*(ustar/uth0-1))
    
    # Grain speed
    
#    ets = np.zeros(s['ustar'].shape)
#    etn = np.zeros(s['ustar'].shape)
#    delhs = np.zeros(s['ustar'].shape)
#    delhn = np.zeros(s['ustar'].shape)
#
#    ix = s['tau'] != 0.
#    ets[ix] = s['taus'][ix]/s['tau'][ix]
#    etn[ix] = s['taun'][ix]/s['tau'][ix]
#    
#
#    
#    delhs[:,:-1] = (s['zb'][:,1:]-s['zb'][:,:-1])/(s['xz'][:,1:]-s['xz'][:,:-1])
#    delhs[:,-1] = delhs[:,-2]
#    delhn[:-1,:] = (s['zb'][1:,:]-s['zb'][:-1,:])/(s['yz'][1:,:]-s['yz'][:-1,:])
#    delhn[-1,:] = delhn[-2,:]
#
#    s['delhs'] = delhs #TEMP
#    s['delhn'] = delhn
#    
#    ets = ets[:,:,np.newaxis].repeat(nf, axis=2)
#    etn = etn[:,:,np.newaxis].repeat(nf, axis=2)
#    delhs = delhs[:,:,np.newaxis].repeat(nf, axis=2)
#    delhn = delhn[:,:,np.newaxis].repeat(nf, axis=2)
#    
#    delhn[:,:,:] = 0.
#    
#    As = np.zeros(ustar.shape)
#    An = np.zeros(ustar.shape)
#    
#    As = np.abs(ets+2*alfa*delhs) #????
#    An = np.abs(etn+2*alfa*delhn)
#    
#    s['As'] = As #TEMP
#    s['An'] = An
#    
#    s['ets'] = ets
#    s['etn'] = etn
#
#    s['uus'][:,:,:] = 0
#    s['uun'][:,:,:] = 0
#    
#    ix = As != 0.
#    s['uus'][ix] = (ueff[ix] - uf/(np.sqrt(2*alfa)*As[ix]))*ets[ix]-np.sqrt(2*alfa)*uf*delhs[ix]/As[ix]
#    
#    ix = etn != 0.
#    s['uun'][ix] = (ueff[ix] - uf/(np.sqrt(2*alfa)*An[ix]))*etn[ix]#-np.sqrt(2*alfa)*uf*delhn[ix]/An[ix]
#    
#    s['uu'] = np.sqrt(s['uus']**2+s['uun']**2)
    
    s['uu'] = (ueff - uf/(np.sqrt(2*alfa)))
    
    ix = ustar != 0.
    
    s['uus'][ix] = s['uu'][ix]*ustars[ix]/ustar[ix]
    s['uun'][ix] = s['uu'][ix]*ustarn[ix]/ustar[ix]
    
#    ix = s['uw'] != 0.
#    
#    s['uus'][ix] = s['uu'][ix]*s['uws'][ix]/s['uw'][ix]
#    s['uun'][ix] = s['uu'][ix]*s['uwn'][ix]/s['uw'][ix]
    
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
        uu = s['uu'][:,:,np.newaxis].repeat(nf, axis=2)


        ix = ustar != 0.
        
        # equilibrium concentration
        
        s['Cu'] = np.zeros(ustar.shape)
        
        if p['method_transport'].lower() == 'bagnold':
            s['Cu'][ix] = np.maximum(0., p['Cb'] * p['rhoa'] / p['g'] \
                                     * (ustar[ix] - s['uth'][ix])**3 / uu[ix])
        elif p['method_transport'].lower() == 'kawamura':
            s['Cu'][ix] = np.maximum(0., p['Cb'] * p['rhoa'] / p['g'] \
                                     * (ustar[ix] + s['uth'][ix])**2 * (ustar[ix] - s['uth'][ix]) / uu[ix])
        elif p['method_transport'].lower() == 'lettau':
            s['Cu'][ix] = np.maximum(0., p['Cb'] * p['rhoa'] / p['g'] \
                                     * (ustar[ix] - s['uth'][ix]) * ustar[ix]**2 / uu[ix])
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

def saturation_factor(s,p):
    
    # THIS PART HAS TO BE CHECKED (BART) - frac
    
#    alpha = 0.35
    gamma = 0.2
    g = 9.81
    
#    frac = 0.
##    
#    Ts_min = 0.0000001
##    Ts_max = 25.
#    
#    nf = p['nfractions']
#    
#    mass = np.zeros(s['uu'].shape)
#    
#    for i in range(nf): #loop over fractions
#        mass[:,:,i] = s['mass'][:,:,0,i]
#        
#        if np.any(s['uth0'][:,:,i]>s['ustar'][:,:]):
#            mass[:,:,i] = 0.
#            
#    for i in range(nf):
#        if np.any(np.max(mass[:,:,:]) == mass[:,:,i]):
#            frac = i
#
#    s['uufrac'][:,:] = s['uu'][:,:,frac]
#    s['uusfrac'][:,:] = s['uus'][:,:,frac]
#    s['uunfrac'][:,:] = s['uun'][:,:,frac]
    
#    uth = s['uth0'][:,:,frac]
    
#    s['T0'] = 2.*alpha*s['ustar']/(g*gamma)
#    s['Ts'] = s['T0']/((s['ustar']/uth)**2.-1.)
#    s['Ts'] = np.maximum(s['Ts'],Ts_min)
    
    s['Ts'] = 1.75*s['ustar']/(g*gamma) #M. Sorensen, Acta Mech., Suppl. 1, 67 1991 FIND PAPER
    s['Ts'] = np.maximum(s['Ts'],0.01)
    
    return s

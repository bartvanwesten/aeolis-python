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

import numpy as np
import logging
import operator

# package modules
import aeolis.shear
from aeolis.utils import *


# initialize logger
logger = logging.getLogger(__name__)


def initialize(s, p):
    '''Initialize wind model

    '''

    # apply wind direction convention
    if isarray(p['wind_file']):
        if p['wind_convention'] == 'cartesian':
            pass
        elif p['wind_convention'] == 'nautical': #CHECK DIRECTIONS!
            p['wind_file'][:,2] = 270.0 - p['wind_file'][:,2]
#            p['wind_file'][:,2] = - p['wind_file'][:,2]
        else:
            logger.log_and_raise('Unknown convention: %s' % p['wind_convention'], exc=ValueError)

    # initialize wind shear model
    if p['process_shear']:
        s['shear'] = aeolis.shear.WindShear(s['x'], s['y'], s['zb'],
                                            L=100., l=10., z0=0.001,
                                            buffer_width=10.)
        
    return s

def interpolate(s, p, t):
    '''Interpolate wind velocity and direction to current time step

    Interpolates the wind time series for velocity and direction to
    the current time step. The cosine and sine of the direction angle
    are interpolated separately to prevent zero-crossing errors. The
    wind velocity is decomposed in two grid components based on the
    orientation of each individual grid cell. In case of a
    one-dimensional model only a single positive component is used.

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters
    t : float
        Current time

    Returns
    -------
    dict
        Spatial grids

    '''
        
    if p['process_wind'] and p['wind_file'] is not None:

        uw_t = p['wind_file'][:,0]
        uw_s = p['wind_file'][:,1]
        uw_d = p['wind_file'][:,2] / 180. * np.pi

        s['uw'][:,:] = interp_circular(t, uw_t, uw_s)
        s['udir'][:,:] = np.arctan2(np.interp(t, uw_t, np.sin(uw_d)),
                                    np.interp(t, uw_t, np.cos(uw_d))) * 180. / np.pi

    s['uws'] = s['uw'] * np.cos(s['udir'] / 180. * np.pi + (s['alfaz'] - 0.5 * np.pi) )
    s['uwn'] = s['uw'] * np.sin(s['udir'] / 180. * np.pi + (s['alfaz'] - 0.5 * np.pi) )
    
    s['uwx'] = s['uw'] * np.cos(s['udir'] / 180. * np.pi)
    s['uwy'] = s['uw'] * np.sin(s['udir'] / 180. * np.pi)

    if p['ny'] == 0:
        s['uwn'][:,:] = 0.
        
    s['uw'] = np.abs(s['uw'])

    # compute saltation velocity
    z = p['z']
    z0 = p['k'] #p['grain_size'][0]/20. #TEMP
    
    # Determine shear velocity
    
    s['ustars'] = s['uws']*p['karman']/np.log(z/z0)
    s['ustarn'] = s['uwn']*p['karman']/np.log(z/z0)
    
    s['ustar'] = np.hypot(s['ustars'],s['ustarn'])
    s['ustar0'] = s['ustar']
    
    # Determine shear stress
    
    s['tau'] = p['rhoa']*s['ustar']**2
    s['tau0'] = s['tau']
    
    s['taus'] = s['tau0']*s['ustars']/s['ustar']
    s['taun'] = s['tau0']*s['ustarn']/s['ustar']
    
#    s['taus0'] = s['taus'].copy
#    s['taun0'] = s['taun'].copy
    
    return s
    
def shear(s,p):
    
    if 'shear' in s.keys() and p['process_shear']:
        
        s['zshear'] = s['zb'].copy()
        
        ix = s['zsep'] > s['zb']
        s['zshear'][ix] = s['zsep'][ix]
        
        # zshear
        
        s['shear'].set_topo(s['zb'])
        s['shear'].set_shear(s['taus'], s['taun'])
        
        s['shear'](u0=s['uw'][0,0],
                   udir=s['udir'][0,0])
        
        s['taus'], s['taun'] = s['shear'].get_shear()
        
#        for j in range(0,p['ny']):
#            avg = np.average(s['dtaus'][j,:2])
#            s['dtaus'][j,:] -= avg
            
#        s['taus'], s['taun'] = s['shear'].add_shear(s['taus'], s['taun'])
        
        # set minimum of taus to zero
        s['tau'] = np.hypot(s['taus'], s['taun'])
        
        ix = s['tau'] != 0.
        
        s['ustar'] = np.sqrt(s['tau'] /p['rhoa'])
        s['ustars'][ix]  = s['ustar'][ix] *s['taus'][ix] /s['tau'][ix] 
        s['ustarn'][ix]  = s['ustar'][ix] *s['taun'][ix] /s['tau'][ix] 
        
        s['zsep'] = s['shear'].get_separation()
#        s['dzsep'] = np.maximum(s['zsep'] - s['zb'], 0.)
        
        # set boundaries
#        s['tau'][:,0] = s['tau0'][:,0]
#        s['taus'][:,0] = s['tau0'][:,0]
#        s['taus'][:,-1] = s['tau0'][:,-1]
#        s['taun'][:,0] = 0.
#        s['taun'][:,-1] = 0.

        # set boundaries (coast - inactive zone)
#        s['taus'][:,:50] = s['tau0'][:,:50]
#        s['tau'][:,:50] = s['tau0'][:,:50]
#        s['taun'][:,:50] = 0.

    return s
#
#def separation(s,p):
#
#    #Calculate delta for relation separation bubble and shear stress
#    
#    x = s['x']
#    dx = x[0,1]-x[0,0]
#    nf = p['nfractions']
##    
#    # ACCORDING CDM
##    m_tau_sepbub = .05
##    slope = np.tan(np.deg2rad(p['Mcr_dyn']))*dx
##    delta = 1./(slope*m_tau_sepbub)
##    
##    s['zsepdelta'] = np.minimum(np.maximum(1. - delta * s['dzsep'], 0.),1.)
#    
##    if p['process_separation']:
##        
##        s['tau'] *= s['zsepdelta']
##    
##        s['taus'] *= s['tau']*ets
##        s['taun'] = s['tau']*etn
#        
#        # uu
#        
##        zsepdelta = np.repeat(s['zsepdelta'][:,:,np.newaxis], nf, axis = 0)
##        s['uu'] *= zsepdelta
##        
##        ets = np.zeros(s['uu'].shape)
##        etn = np.zeros(s['uu'].shape)
##            
##        ix = s['uu'] != 0.
##        ets[ix] = s['uus'][ix]/s['uu'][ix]
##        etn[ix] = s['uun'][ix]/s['uu'][ix]
##        
##        s['uus'] = s['uu']*ets
##        s['uun'] = s['uu']*etn
#        
#    s['ustar'] = np.sqrt(s['tau']/p['rhoa'])
#    s['ustars'] = s['ustar']*s['taus']/s['tau']
#    s['ustarn'] = s['ustar']*s['taun']/s['tau']
#    
#    return s

def get_velocity_at_height(u, z, z0, z1=None):
    '''Compute shear velocity from wind velocity following Prandl-Karman's Law of the Wall

    Parameters
    ----------
    u : numpy.ndarray
        Spatial wind field
    z : float
        Height above bed where ``u`` is measured
    z0 : float
        Roughness length
    z1 : float, optional
        Height above bed for which to return wind speeds.
        Returns wind shear if not given.

    Returns
    -------
    numpy.ndarray
        Array of size ``u`` with wind speeds at height ``z1``

    '''

    tau = .41 / np.log(z / z0) * u

    if z1 is None:
        return tau
    else:
        return tau * np.log(z1 / z0) / .41
    
def filter_low(s, p, par, direction, Cut):
    
    nx = p['nx'] + 1
    ny = p['ny'] + 1
    
    parfft = np.fft.fft2(s[par])
    dk = 2.0 * np.pi / np.max(s[direction])
    fac = np.exp(-((dk*s[direction])**2.)/(2.*Cut**2.))
    parfft *= fac
    s[par] = np.real(np.fft.ifft2(parfft))
    
    if par == 'tau' or par == 'taunosep' or par == 'dzsep':
        
        zerostate   = np.zeros((ny, nx))
        avg         = np.zeros((ny, nx))
        zeroavg     = np.zeros((ny, nx))
        
        if par == 'dzsep':
            zerostate[:,:] = 0. # adjust this
    
        for j in range(0,p['ny']):
            avg = np.average(s[par][j,:10])
            zeroavg = np.average(zerostate[j,:])
            s[par][j,:]+=(zeroavg-avg)
        for i in range(0,p['nx']):
            avg = np.average(s[par][:10,i])
            zeroavg = np.average(zerostate[:,i])
            s[par][:,i]+=(zeroavg-avg)
    
    return s

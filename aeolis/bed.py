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
    '''Initialize bathymetry and bed composition

    Initialized bathymetry, computes cell sizes and orientation, bed
    layer thickness and bed composition.

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
    
    # get model dimensions # Method suitable for non-equidistant grid, not for curvilinear!
    ny = p['ny']
    nl = p['nlayers']
    nf = p['nfractions']

    # initialize x-dimension
    s['x'][:,:] = p['xgrid_file']
    
    # initialize distance between cell centers 
    s['dsc'][:-1,1:-1] = np.diff(s['x'], axis=1)
    s['dsc'][:-1,0] = s['dsc'][:-1,1]
    s['dsc'][:-1,-1] = s['dsc'][:-1,-2]
    s['dsc'][-1,:] = s['dsc'][-2,:]
    
    # initialize corners xc
    s['xc'][1:,:-1] = s['x']+s['dsc'][1:,:-1]*0.5
    s['xc'][0,:-1] = s['x'][0,:]-s['dsc'][0,:-1]*0.5
    s['xc'][:,-1] = s['xc'][:,-1]
    
    # initialize cell dimensions
    s['ds'][:,:] = np.diff(s['xc'][:-1,:], axis=1)
    
    # initialize y-dimension
    if ny == 0:
        s['y'][:,:] = 0.
        s['dn'][:,:] = 1.
        s['dnc'][:,:] = 1.
        s['alfa'][:,:] = 0.
    else:
        # initialize y-dimension
        s['y'][:,:] = p['ygrid_file']
        
        # initialize distance between cell centers 
        s['dnc'][1:-1,:-1] = np.diff(s['y'], axis=0)
        s['dnc'][0,:-1] = s['dnc'][1,:-1]
        s['dnc'][-1,:-1] = s['dnc'][-2,:-1]
        s['dnc'][:,-1] = s['dnc'][:,-2]
        
        # initialize corners yc 
        s['yc'][:-1,1:] = s['y']+s['dnc'][:-1,1:]*0.5
        s['yc'][:-1,0] = s['y'][:,0]-s['dnc'][:-1,0]*0.5
        s['yc'][-1,:] = s['yc'][-1,:]
        
        # initialize cell dimensions
        s['dn'][:,:] = np.diff(s['yc'][:,:-1], axis=0)

        s['alfa'][1:-1,:] = np.arctan2(s['x'][2:,:] - s['x'][:-2,:],
                                       s['y'][2:,:] - s['y'][:-2,:])
        s['alfa'][0,:] = s['alfa'][1,:]
        s['alfa'][-1,:] = s['alfa'][-2,:]

    # compute cell areas
    s['dsdn'][:,:] = s['ds'] * s['dn']
    s['dsdni'][:,:] = 1. / s['dsdn']
    
    s['dsdnc'][:,:] = s['dsc'] * s['dnc']
    s['dsdnci'][:,:] = 1. / s['dsdnc']

    # initialize bathymetry
    s['zb'][:,:] = p['bed_file']

    # initialize bed layers
    s['thlyr'][:,:,:] = p['layer_thickness']

    # initialize bed composition
    if p['bedcomp_file'] is None:
        gs = makeiterable(p['grain_dist'])
        gs = gs / np.sum(gs)
        for i in range(nl):
            for j in range(nf):
                s['mass'][:,:,i,j] = p['rhop'] * (1. - p['porosity']) \
                                     * s['thlyr'][:,:,i] * gs[j]
    else:
        s['mass'][:,:,:,:] = p['bedcomp_file'].reshape(s['mass'].shape)                

    # initialize masks
    for k, v in p.items():
        if k.endswith('_mask'):
            if v is None:
                s[k] = 1.
            else:
                s[k] = v.reshape(s['zb'].shape)

    # initialize threshold
    if p['threshold_file'] is not None:
        s['uth'] = p['threshold_file'][:,:,np.newaxis].repeat(nf, axis=-1)
        
    return s


def update(s, p):
    '''Update bathymetry and bed composition

    Update bed composition by moving sediment fractions between bed
    layers. The total mass in a single bed layer does not change as
    sediment removed from a layer is repleted with sediment from
    underlying layers. Similarly, excess sediment added in a layer is
    moved to underlying layers in order to keep the layer mass
    constant. The lowest bed layer exchanges sediment with an infinite
    sediment source that follows the original grain size distribution
    as defined in the model configuration file by ``grain_size`` and
    ``grain_dist``. The bathymetry is updated following the
    cummulative erosion/deposition over the fractions if ``bedupdate``
    is ``True``.

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

    nx = p['nx']
    ny = p['ny']
    nl = p['nlayers']
    nf = p['nfractions']

    # determine net erosion
    pickup = s['pickup'].reshape((-1,nf))

    # determine total mass that should be exchanged between layers
    dm = -np.sum(pickup, axis=-1, keepdims=True).repeat(nf, axis=-1)
    
    # get erosion and deposition cells
    ix_ero = dm[:,0] < 0.
    ix_dep = dm[:,0] > 0.
    
    # reshape mass matrix
    m = s['mass'].reshape((-1,nl,nf))

    # negative mass may occur in case of deposition due to numerics,
    # which should be prevented
    m, dm, pickup = prevent_negative_mass(m, dm, pickup)
    
    # determine weighing factors
    d = normalize(m, axis=2)
    
    # move mass among layers
    m[:,0,:] -= pickup
    for i in range(1,nl):
        m[ix_ero,i-1,:] -= dm[ix_ero,:] * d[ix_ero,i,:]
        m[ix_ero,i,  :] += dm[ix_ero,:] * d[ix_ero,i,:]
        m[ix_dep,i-1,:] -= dm[ix_dep,:] * d[ix_dep,i-1,:]
        m[ix_dep,i,  :] += dm[ix_dep,:] * d[ix_dep,i-1,:]
    m[ix_dep,-1,:] -= dm[ix_dep,:] * d[ix_dep,-1,:]
    m[ix_ero,-1,:] -= dm[ix_ero,:] * normalize(p['grain_dist'])[np.newaxis,:].repeat(np.sum(ix_ero), axis=0)

    # remove tiny negatives
    m = prevent_tiny_negatives(m, p['max_error'])

    # warn if not all negatives are gone
    if m.min() < 0:
        logger.warning(format_log('Negative mass',
                                  nrcells=np.sum(np.any(m<0., axis=-1)),
                                  minvalue=m.min(),
                                  minwind=s['uw'].min(),
                                  time=p['_time']))
        
    # reshape mass matrix
    s['mass'] = m.reshape((ny+1,nx+1,nl,nf))

    # update bathy
    if p['process_bedupdate']:
        dz = dm[:,0].reshape((ny+1,nx+1)) / (p['rhop'] * (1. - p['porosity']))
        s['zb'] += dz
        s['zs'] += dz

    return s


def prevent_negative_mass(m, dm, pickup):
    '''Handle situations in which negative mass may occur due to numerics

    Negative mass may occur by moving sediment to lower layers down to
    accomodate deposition of sediments. In particular two cases are
    important:

    #. A net deposition cell has some erosional fractions.

       In this case the top layer mass is reduced according to the
       existing sediment distribution in the layer to accomodate
       deposition of fresh sediment. If the erosional fraction is
       subtracted afterwards, negative values may occur. Therefore the
       erosional fractions are subtracted from the top layer
       beforehand in this function. An equal mass of deposition
       fractions is added to the top layer in order to keep the total
       layer mass constant. Subsequently, the distribution of the
       sediment to be moved to lower layers is determined and the
       remaining deposits are accomodated.

    #. Deposition is larger than the total mass in a layer.

       In this case a non-uniform distribution in the bed may also
       lead to negative values as the abundant fractions are reduced
       disproportionally as sediment is moved to lower layers to
       accomodate the deposits. This function fills the top layers
       entirely with fresh deposits and moves the existing sediment
       down such that the remaining deposits have a total mass less
       than the total bed layer mass. Only the remaining deposits are
       fed to the routine that moves sediment through the layers.

    Parameters
    ----------
    m : np.ndarray
        Sediment mass in bed (nx*ny, nl, nf)
    dm : np.ndarray
        Total sediment mass exchanged between layers (nx*ny, nf)
    pickup : np.ndarray
        Sediment pickup (nx*ny, nf)

    Returns
    -------
    np.ndarray
        Sediment mass in bed (nx*ny, nl, nf)
    np.ndarray
        Total sediment mass exchanged between layers (nx*ny, nf)
    np.ndarray
        Sediment pickup (nx*ny, nf)

    Note
    ----
    The situations handled in this function can also be prevented by
    reducing the time step, increasing the layer mass or increasing
    the adaptation time scale.

    '''

    nl = m.shape[1]
    nf = m.shape[2]

    ###
    ### case #1: deposition cells with some erosional fractions
    ###
    
    ix_dep = dm[:,0] > 0.
    
    # determine erosion and deposition fractions per cell
    ero =  np.maximum(0., pickup)
    dep = -np.minimum(0., pickup)

    # determine gross erosion
    erog = np.sum(ero, axis=1, keepdims=True).repeat(nf, axis=1)

    # determine net deposition cells with some erosional fractions
    ix = ix_dep & (erog[:,0] > 0)

    # remove erosional fractions from pickup and remove an equal mass
    # of accretive fractions from the pickup, adapt sediment exchange
    # mass and bed composition accordingly
    if np.any(ix):
        d = normalize(dep, axis=1)
        ddep = erog[ix,:] * d[ix,:]
        pickup[ix,:] = -dep[ix,:] + ddep
        dm[ix,:] = -np.sum(pickup[ix,:], axis=-1, keepdims=True).repeat(nf, axis=-1)
        m[ix,0,:] -= ero[ix,:] - ddep # FIXME: do not use deposition in normalization

    ###
    ### case #2: deposition cells with deposition larger than the mass present in the top layer
    ###

    mx = m[:,0,:].sum(axis=-1, keepdims=True)

    # determine deposition in terms of layer mass (round down)
    n = dm[:,:1] // mx

    # determine if deposition is larger than a sinle layer mass
    if np.any(n > 0):

        # determine distribution of deposition
        d = normalize(pickup, axis=1)

        # walk through layers from top to bottom
        for i in range(nl):

            ix = (n > i).flatten()
            if not np.any(ix):
                break

            # move all sediment below current layer down one layer
            m[ix,(i+1):,:] = m[ix,i:-1,:]

            # fill current layer with deposited sediment
            m[ix,i,:] = mx[ix,:].repeat(nf, axis=1) * d[ix,:]

            # remove deposited sediment from pickup
            pickup[ix,:] -= m[ix,i,:]

        # discard any remaining deposits at locations where all layers
        # are filled with fresh deposits
        ix = (dm[:,:1] > mx).flatten()
        if np.any(ix):
            pickup[ix,:] = 0.

        # recompute sediment exchange mass
        dm[ix,:] = -np.sum(pickup[ix,:], axis=-1, keepdims=True).repeat(nf, axis=-1)

    return m, dm, pickup


def mixtoplayer(s, p):
    '''Mix grain size distribution in top layer of the bed

    Simulates mixing of the top layers of the bed by wave action. The
    wave action is represented by a local wave height maximized by a
    maximum wave hieght over depth ratio ``gamma``. The mixing depth
    is a fraction of the local wave height indicated by
    ``facDOD``. The mixing depth is used to compute the number of bed
    layers that should be included in the mixing. The grain size
    distribution in these layers is then replaced by the average grain
    size distribution over these layers.

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

    if p['process_mixtoplayer']:
        
        # get model dimensions
        nx = p['nx']+1
        ny = p['ny']+1
        nl = p['nlayers']
        nf = p['nfractions']

        # compute depth of disturbence for each cell and repeat for each layer
        DOD = p['facDOD'] * s['Hs']

        # compute ratio total layer thickness and depth of disturbance 
        ix = DOD > 0.
        f = np.ones(DOD.shape)
        f[ix] = np.minimum(1., s['thlyr'].sum(axis=2)[ix] / DOD[ix])

        # correct shapes
        DOD = DOD[:,:,np.newaxis].repeat(nl, axis=2)
        f = f[:,:,np.newaxis].repeat(nl, axis=2)

        # determine what layers are above the depth of disturbance
        ix = (s['thlyr'].cumsum(axis=2) <= DOD) & (DOD > 0.)
        ix = ix[:,:,:,np.newaxis].repeat(nf, axis=3)
        f = f[:,:,:,np.newaxis].repeat(nf, axis=3)
        
        # average mass over layers
        if np.any(ix):
            ix[:,:,0,:] = True # at least mix the top layer
            mass = s['mass'].copy()
            mass[~ix] = np.nan
            
            gd = normalize(p['grain_dist']) * p['rhop'] * (1. - p['porosity'])
            gd = gd.reshape((1,1,1,-1)).repeat(ny, axis=0) \
                                       .repeat(nx, axis=1) \
                                       .repeat(nl, axis=2)

            mass1 = np.nanmean(mass, axis=2, keepdims=True).repeat(nl, axis=2)
            mass2 = gd * s['thlyr'][:,:,:,np.newaxis].repeat(nf, axis=-1)
            mass = mass1 * f + mass2 * (1. - f)
        
            s['mass'][ix] = mass[ix]
            
    return s

# NEW!
def avalanche(s, p):
    
    if p['process_avalanche']:
    
        nx = p['nx']+1
        ny = p['ny']+1
        
        x = s['x']
        y = s['y']
        
        dsc = s['dsc']
        dnc = s['dnc']
        
        # Calculation of ratio dsdn
        
        dsdn = np.repeat(s['dsdn'][:,:,np.newaxis], 8, axis = 2)
        
        dsdn_neighbour = np.zeros((ny,nx,8))
        dsdn_neighbour[1:-1,:-2,0] = dsdn[1:-1,1:-1,0]
        dsdn_neighbour[1:-1,2:,1] = dsdn[1:-1,1:-1,0]
        dsdn_neighbour[:-2,1:-1,2] = dsdn[1:-1,1:-1,0]
        dsdn_neighbour[2:,1:-1,3] = dsdn[1:-1,1:-1,0]
        dsdn_neighbour[2:,:-2,4] = dsdn[1:-1,1:-1,0]
        dsdn_neighbour[2:,2:,5] = dsdn[1:-1,1:-1,0]
        dsdn_neighbour[:-2,2:,6] = dsdn[1:-1,1:-1,0]
        dsdn_neighbour[:-2,:-2,7] = dsdn[1:-1,1:-1,0]
        
        dsdn_fac = np.divide(dsdn, dsdn_neighbour, out=np.zeros_like(dsdn), where = dsdn_neighbour != 0)
        
        #parameters
        Mcr_stat = 34
        Mcr_dyn = 33
        
        # Calculation of dh
        dh_stat =        np.zeros((ny,nx,8))+1000 # 1000, because np.inf causes errors later on
        dh_dyn =         np.zeros((ny,nx,8))+1000
        #statdyndiff =    np.zeros((NY+1,NX+1,8))+1000
        
        # Calculation of dh_stat
        dh_stat[:,:-1,0] = np.tan(Mcr_stat*(np.pi/180.))*(dsc[:-1,1:-1]) #Negative X-direction
        dh_stat[:,1:,1]  = np.tan(Mcr_stat*(np.pi/180.))*(dsc[:-1,1:-1]) #Positive X-direction
        dh_stat[:-1,:,2] = np.tan(Mcr_stat*(np.pi/180.))*(dnc[1:-1,:-1]) #Negative Y-direction
        dh_stat[1:,:,3]  = np.tan(Mcr_stat*(np.pi/180.))*(dnc[1:-1,:-1]) #Positive Y-direction
        dh_stat[1:,:-1,4]  = np.tan(Mcr_stat*(np.pi/180.))*((dsc[1:-1,1:-1])**2.+(dnc[1:-1,1:-1])**2.)**0.5 #Negative X-direction and Positive Y-direction
        dh_stat[1:,1:,5]   = np.tan(Mcr_stat*(np.pi/180.))*((dsc[1:-1,1:-1])**2.+(dnc[1:-1,1:-1])**2.)**0.5 #Positive X-direction and Positive Y-direction
        dh_stat[:-1,1:,6]  = np.tan(Mcr_stat*(np.pi/180.))*((dsc[1:-1,1:-1])**2.+(dnc[1:-1,1:-1])**2.)**0.5 #Positive X-direction and Negative Y-direction
        dh_stat[:-1,:-1,7] = np.tan(Mcr_stat*(np.pi/180.))*((dsc[1:-1,1:-1])**2.+(dnc[1:-1,1:-1])**2.)**0.5 #Negative X-direction and Negative Y-direction
        
        # Calculation of dh_stat
        dh_dyn[:,:-1,0] = np.tan(Mcr_dyn*(np.pi/180.))*(dsc[:-1,1:-1]) #Negative X-direction
        dh_dyn[:,1:,1]  = np.tan(Mcr_dyn*(np.pi/180.))*(dsc[:-1,1:-1]) #Positive X-direction
        dh_dyn[:-1,:,2] = np.tan(Mcr_dyn*(np.pi/180.))*(dnc[1:-1,:-1]) #Negative Y-direction
        dh_dyn[1:,:,3]  = np.tan(Mcr_dyn*(np.pi/180.))*(dnc[1:-1,:-1]) #Positive Y-direction
        dh_dyn[1:,:-1,4]  = np.tan(Mcr_dyn*(np.pi/180.))*((dsc[1:-1,1:-1])**2.+(dnc[1:-1,1:-1])**2.)**0.5 #Negative X-direction and Positive Y-direction
        dh_dyn[1:,1:,5]   = np.tan(Mcr_dyn*(np.pi/180.))*((dsc[1:-1,1:-1])**2.+(dnc[1:-1,1:-1])**2.)**0.5 #Positive X-direction and Positive Y-direction
        dh_dyn[:-1,1:,6]  = np.tan(Mcr_dyn*(np.pi/180.))*((dsc[1:-1,1:-1])**2.+(dnc[1:-1,1:-1])**2.)**0.5 #Positive X-direction and Negative Y-direction
        dh_dyn[:-1,:-1,7] = np.tan(Mcr_dyn*(np.pi/180.))*((dsc[1:-1,1:-1])**2.+(dnc[1:-1,1:-1])**2.)**0.5 #Negative X-direction and Negative Y-direction
        
        # Calculation of statdyndiff
        statdyndiff = dh_stat - dh_dyn
        
        # Reset counter for while-loop
        count=-1
        
        # Initialize different bathymetries
        zb = s['zb']
        zb_center = np.zeros((ny,nx,8))
        zb_neighbour = np.zeros((ny,nx,8))
        
        # Start of the while-loop
        for i in range(p['max_iter']):
                    
            count+=1
            
            #ZB neighbour
            zb_center = np.repeat(zb[:,:,np.newaxis], 8, axis = 2)
            #Negative X-direction
            zb_neighbour[:,:-1,0] = zb[:,1:]
            zb_neighbour[:,-1,0] = zb_neighbour[:,-2,0]
            #Positive X-direction
            zb_neighbour[:,1:,1] = zb[:,:-1]
            zb_neighbour[:,0,1] = zb_neighbour[:,1,1]
            #Negative Y-direction
            zb_neighbour[:-1,:,2] = zb[1:,:]
            zb_neighbour[-1,:,2] = zb_neighbour[-2,:,2]
            #Positive Y-direction
            zb_neighbour[1:,:,3] = zb[:-1,:]
            zb_neighbour[0,:,3] = zb_neighbour[1,:,3]
            #XY
            zb_neighbour[1:,:-1,4] = zb[:-1,1:]
            zb_neighbour[0,:,4] = zb_neighbour[1,:,4]
            zb_neighbour[:,-1,4] = zb_neighbour[:,-2,4]
            #XY
            zb_neighbour[1:,1:,5] = zb[:-1,:-1]
            zb_neighbour[0,:,5] = zb_neighbour[1,:,5]
            zb_neighbour[:,0,5] = zb_neighbour[:,1,5]
            #XY
            zb_neighbour[:-1,1:,6] = zb[1:,:-1]
            zb_neighbour[-1,:,6] = zb_neighbour[-2,:,6]
            zb_neighbour[:,0,6] = zb_neighbour[:,1,6]
            #XY
            zb_neighbour[:-1,:-1,7] = zb[1:,1:]
            zb_neighbour[-1,:,7] = zb_neighbour[-2,:,7]
            zb_neighbour[:,-1,7] = zb_neighbour[:,-2,7]    
            
            # Reset arrays
            total   = np.zeros((ny, nx))
            surplus = np.zeros((ny, nx, 8))
            flux    = np.zeros((ny, nx, 8))
            
            # Calculate the surpluses, filter negative values and add the statdyndiff
            surplus = (zb_center - zb_neighbour) - dh_stat
            surplus = surplus.clip(0.)
            surplus += (surplus>0) * statdyndiff
            
            # Sum all the surplusses
            total = surplus.sum(axis=2)
            total = np.repeat(total[:,:,np.newaxis], 8, axis = 2)
            
            if np.sum(total) == 0:
                break
            else:
            
                flux = surplus * np.divide(surplus, total + surplus, out=np.zeros_like(surplus), where = total + surplus != 0)
                totalflux = flux.sum(axis=2)
                maxflux = flux.max(axis=2)
                maxsurplus = surplus.max(axis=2)
                
                #Fluxnorm   
                a1 = maxflux
                b1 = totalflux
                c1 = np.divide(a1, b1, out=np.zeros_like(a1), where = b1 != 0)
                a2 = maxsurplus
                b2 = 1 + c1
                c2 = np.divide(a2, b2, out=np.zeros_like(a2), where = b2 != 0)
                a3 = c2
                b3 = totalflux
                fluxnorm = np.divide(a3, b3, out=np.zeros_like(a3), where = b3 != 0)
                
                #restribute sediment
                zb          -= fluxnorm * totalflux# * dsdn_fac
                zb[:,:-1]   += fluxnorm[:,1:] * flux[:,1:,1] * dsdn_fac[:,1:,1]
                zb[:,1:]    += fluxnorm[:,:-1] * flux[:,:-1,0] * dsdn_fac[:,:-1,0]
                zb[:-1,:]   += fluxnorm[1:,:] * flux[1:,:,3] * dsdn_fac[1:,:,3]
                zb[1:,:]    += fluxnorm[:-1,:] * flux[:-1,:,2] * dsdn_fac[:-1,:,2]
                zb[1:,:-1]  += fluxnorm[:-1,1:] * flux[:-1,1:,6] * dsdn_fac[:-1,1:,6]
                zb[1:,1:]   += fluxnorm[:-1,:-1] * flux[:-1,:-1,7] * dsdn_fac[:-1,:-1,7]
                zb[:-1,1:]  += fluxnorm[1:,:-1] * flux[1:,:-1,4] * dsdn_fac[1:,:-1,4]
                zb[:-1,:-1] += fluxnorm[1:,1:] * flux[1:,1:,5] * dsdn_fac[1:,1:,5]
                
                #Boundary
                zb[0,:]  = zb[1,:]
                zb[-1,]  = zb[-2,:]
                zb[:,0]  = zb[:,1]
                zb[:,-1] = zb[:,-2]
        
        s['zb']=zb
    
    return s

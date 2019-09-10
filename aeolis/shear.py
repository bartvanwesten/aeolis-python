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

import logging
import numpy as np
import scipy.special
import scipy.interpolate
import matplotlib.pyplot as plt


# initialize logger
logger = logging.getLogger(__name__)


class WindShear:
    '''Class for computation of 2DH wind shear perturbations over a topography.
        
    The class implements a 2D FFT solution to the wind shear
    perturbation on curvilinear grids.  As the FFT solution is only
    defined on an equidistant rectilinear grid with circular boundary
    conditions that is aligned with the wind direction, a rotating
    computational grid is automatically defined for the computation.
    The computational grid is extended in all directions using a
    logistic sigmoid function as to ensure full coverage of the input
    grid for all wind directions, circular boundaries and preservation
    of the alongshore uniformity.  An extra buffer distance can be
    used as to minimize the disturbence from the borders in the input
    grid.  The results are interpolated back to the input grid when
    necessary.

    Frequencies related to wave lengths smaller than a computational
    grid cell are filtered from the 2D spectrum of the topography
    using a logistic sigmoid tapering. The filtering aims to minimize
    the disturbance as a result of discontinuities in the topography
    that may physically exists, but cannot be solved for in the
    computational grid used.

    Example
    -------
    >>> w = WindShear(x, y, z)
    >>> w(u0=10., udir=30.).add_shear(taux, tauy)

    Notes
    -----
    To do:

    * Actual resulting values are still to be compared with the results
       from Kroy et al. (2002)
    * Grid interpolation can still be optimized
    * Separation bubble is still to be implemented
    * Avalanching is still to be implemented

    '''

    
    igrid = {}
    cgrid = {}
    istransect = False
    
    
    def __init__(self, x, y, z, dx=.5, dy=.5,
                 buffer_width=100., buffer_relaxation=None,
                 L=100., d=0.025, l=10.):
        '''Class initialization
            
        Parameters
        ----------
        x : numpy.ndarray
            2D array with x-coordinates of input grid
        y : numpy.ndarray
            2D array with y-coordinates of input grid
        z : numpy.ndarray
            2D array with topography of input grid
        dx : float, optional
            Grid spacing in x dimension of computational grid
            (default: 1)
        dy : float, optional
            Grid spacing of y dimension of computational grid
            (default: 1)
        buffer_width : float, optional
            Width of buffer distance between input grid boundary and
            computational grid boundary (default: 100)
        buffer_relaxation : float, optional
            Relaxation of topography in buffer from input grid
            boundary to computational grid boundary (default:
            buffer_width / 4)
        L : float, optional
            Length scale of topographic features (default: 100)
        z0 : float, optional
            Aerodynamic roughness (default: .001)
        l : float, optional
            Height of inner layer (default: 10)

        '''
        
        if buffer_relaxation is None:
            buffer_relaxation = buffer_width / 4.

        if z.shape[0] == 1:
            self.istransect = True
        
        self.igrid = dict(x = x,
                          y = y,
                          z = z)

        self.irgrid = dict(x = x,
                           y = y,
                           z = z)
            
        self.cgrid = dict(dx = dx,
                          dy = dy)
                          
        self.buffer_width = buffer_width
        self.buffer_relaxation = buffer_relaxation
                          
        self.L = L
        self.z0 = d[0]/20. # z0 # ADJUST THIS TO MULTIPLE FACTIONS!
        self.l = l
                          
        self.set_computational_grid()


    def __call__(self, u0, udir, process_separation):
        '''Compute wind shear for given wind speed and direction
        
        Parameters
        ----------
        u0 : float
            Free-flow wind speed
        udir : float
            Wind direction in degrees
        
        '''
        u1 = np.zeros(u0.shape)
        
        ix = u0 > 0.
        u1[ix] = u0[ix]/u0[ix] * 0.76 #why??

        udir = 90 + udir

        # Populate computational-grid with topography
        self.populate_computational_grid(udir)

        # Compute separation bubble
        if process_separation:
            zsep = self.separation()
            z_origin = self.cgrid['z'].copy()
            self.cgrid['z'] = np.maximum(self.cgrid['z'], zsep)
        
        # Compute shear stresses on computational grid
        self.compute_shear(u1)
                    
        gc = self.cgrid
        gi = self.igrid
        
        gc['dtaux'] = np.maximum(gc['dtaux'],-1.)
        
        # gc['taux'], gc['tauy'] = self.rotate(gc['taux'], gc['tauy'], udir)

        # Add shear and multiply shear with separation
        self.add_shear()

        if process_separation:
            self.cgrid['dzsep'] = self.cgrid['z'] - z_origin
            self.separation_shear(self.cgrid['dzsep'])
        
        # Rotation
        gi['x'], gi['y'] =          self.rotate(gi['x'],    gi['y'],    -udir, origin=(self.x0, self.y0))
        gc['x'], gc['y'] =          self.rotate(gc['x'],    gc['y'],    -udir, origin=(self.x0, self.y0))
        gc['taux'], gc['tauy'] =    self.rotate(gc['taux'], gc['tauy'], -udir)

        # Interpolate real grid with shear stresses and separation bubble
        self.igrid['taux'] = self.interpolate(gc['x'], gc['y'], gc['taux'], gi['x'], gi['y'])
        self.igrid['tauy'] = self.interpolate(gc['x'], gc['y'], gc['tauy'], gi['x'], gi['y'])
        self.igrid['dzsep'] = self.interpolate(gc['x'], gc['y'], gc['dzsep'], gi['x'], gi['y'])

        # Filled rotational grid, now rotate it back to original direction

        gi['taux'], gi['tauy'] = self.rotate(gi['taux'], gi['tauy'], udir)
        gi['x'], gi['y'] = self.rotate(gi['x'], gi['y'], udir, origin=(self.x0, self.y0))

        # PLOTTING! -------------------------------------------------------
        # d = 10
        # plt.pcolormesh(gi['x'], gi['y'], gi['dzsep'])
        # plt.quiver(gi['x'][::d, ::d], gi['y'][::d, ::d],
        #           gi['taux'][::d, ::d], gi['tauy'][::d, ::d])
        # plt.colorbar()
        # plt.show()
        # -------------------------------------------------

        gc['x'], gc['y'] = self.rotate(gc['x'], gc['y'], udir, origin=(self.x0, self.y0))

        return self


    def get_shear(self):
        '''Returns wind shear perturbation field
        
        Returns
        -------
        dtaux : numpy.ndarray
            Wind shear perturbation in x-direction
        dtauy : numpy.ndarray
            Wind shear perturbation in y-direction
        
        '''

        taux = self.igrid['taux']
        tauy = self.igrid['tauy']
            
        return taux, tauy
    
    def get_separation(self):
        '''Returns wind shear perturbation field
        
        Returns
        -------
        dtaux : numpy.ndarray
            Wind shear perturbation in x-direction
        dtauy : numpy.ndarray
            Wind shear perturbation in y-direction
        
        '''

        zsep = self.igrid['dzsep']
            
        return zsep
    
    def get_hs(self):
        '''Returns wind shear perturbation field
        
        Returns
        -------
        dtaux : numpy.ndarray
            Wind shear perturbation in x-direction
        dtauy : numpy.ndarray
            Wind shear perturbation in y-direction
        
        '''

        hs0 = self.hs0
        hs1 = self.hs1
            
        return hs0, hs1
        
        
    def add_shear(self):
        '''Add wind shear perturbations to a given wind shear field
        
        Parameters
        ----------
        taux : numpy.ndarray
            Wind shear in x-direction
        tauy : numpy.ndarray
            Wind shear in y-direction

        Returns
        -------
        taux : numpy.ndarray
            Wind shear including perturbations in x-direction
        tauy : numpy.ndarray
            Wind shear including perturbations in y-direction
        
        '''

        tau = np.sqrt(self.cgrid['taux']**2 + self.cgrid['tauy']**2)
        ix = tau > 0.

        dtaux = self.cgrid['dtaux']
        dtauy = self.cgrid['dtauy']
        
        self.cgrid['taux'][ix] = tau[ix] * (self.cgrid['taux'][ix] / tau[ix] + dtaux[ix])
        self.cgrid['tauy'][ix] = tau[ix] * (self.cgrid['tauy'][ix] / tau[ix] + dtauy[ix])

        return


    def set_topo(self, z):
        '''Update topography

        Parameters
        ----------
        z : numpy.ndarray
            2D array with topography of input grid

        '''

        self.igrid['z'] = z

        return self
        
    def set_shear(self, taus, taun):
        '''Update shear

        Parameters
        ----------
        z : numpy.ndarray
            2D array with topography of input grid

        '''
        self.igrid['taux'] = taus
        self.igrid['tauy'] = taun

        return self
    
    def populate_computational_grid(self, alpha):
        '''Interpolate input topography to computational grid
            
        Rotates computational grid to current wind direction and
        interpolates the input topography to the rotated grid. Any
        grid cells that are not covered by the input grid are filled
        using a sigmoid function.
            
        Parameters
        ----------
        alpha : float
            Rotation angle in degrees

        '''
        
        gc = self.cgrid
        gi = self.igrid
        
        # Add buffer zone around grid

        dxi = gi['x'][1,1]-gi['x'][0,0]
        dyi = gi['y'][1,1]-gi['y'][0,0]

        buf = 200 # amount of cells

        xi, yi = np.meshgrid(np.linspace(gi['x'][0,0]-buf*dxi,gi['x'][-1,-1]+buf*dxi,gi['x'].shape[1]+2*buf),
                            np.linspace(gi['y'][0,0]-buf*dyi,gi['y'][-1,-1]+buf*dyi,gi['y'].shape[0]+2*buf))

        zi = np.zeros((xi.shape))

        zi[buf:-buf,buf:-buf] = gi['z']

        # Filling bufferzone (edges)

        zi[buf:-buf,:buf] = np.repeat(zi[buf:-buf,buf+1][:,np.newaxis], buf, axis = 1)
        zi[buf:-buf,-buf:] = np.repeat(zi[buf:-buf,-buf-1][:,np.newaxis], buf, axis = 1)

        zi[:buf,buf:-buf] = np.repeat(zi[buf+1,buf:-buf][np.newaxis], buf, axis = 0)
        zi[-buf:,buf:-buf] = np.repeat(zi[-buf-1,buf:-buf][np.newaxis], buf, axis = 0)

        # Filling bufferzone (corners)

        zi[:buf,:buf] = zi[buf+1,buf+1]
        zi[-buf:,:buf] = zi[-buf-1,buf+1]
        zi[:buf,-buf:] = zi[buf+1,-buf-1]
        zi[-buf:,-buf:] = zi[-buf-1,-buf-1]

        # Rotate computational

        xc, yc = self.rotate(gc['xi'], gc['yi'], alpha, origin=(self.x0, self.y0))
        zc = self.interpolate(xi, yi, zi, xc, yc)

        tauxc = self.interpolate(gi['x'], gi['y'], gi['taux'], xc, yc)
        tauyc = self.interpolate(gi['x'], gi['y'], gi['tauy'], xc, yc)

        self.cgrid['z'] = zc
        self.cgrid['taux'] = tauxc
        self.cgrid['tauy'] = tauyc
        self.cgrid['x'] = xc
        self.cgrid['y'] = yc
        
        # NO INTERPOLATION
        #
        # self.cgrid['z'] = gi['z']
        # self.cgrid['taux'] = gi['taux']
        # self.cgrid['tauy'] = gi['tauy']
        # self.cgrid['x'] = gi['x']
        # self.cgrid['y'] = gi['y']
        
    def compute_shear(self, u0):
        '''Compute wind shear perturbation for given free-flow wind speed on computational grid
        
        Parameters
        ----------
        u0 : float
            Free-flow wind speed
        nfilter : 2-tuple
            Wavenumber range used for logistic sigmoid filter. See
            :func:`filter_highfrequencies`

        '''
            
        kappa = 0.41 #self.p['karman']
        
        g = self.cgrid
                
        if u0 == 0.:
            self.cgrid['dtaux'] = np.zeros(g['z'].shape)
            self.cgrid['dtauy'] = np.zeros(g['z'].shape)
            return
                                
        ny, nx = g['z'].shape
        kx, ky = np.meshgrid(2. * np.pi * np.fft.fftfreq(nx+1, g['dx'])[1:],
                             2. * np.pi * np.fft.fftfreq(ny+1, g['dy'])[1:])
        
        hs = np.fft.fft2(g['z'])
        
        # Filter
        hs = self.filter_highfrequenies(kx, ky, hs, (1.5, 6), 0.01)

        # 1. Auxiliary variables
        #-----------------------
        # 1.1 Mean lengthscale
        L = np.sum(np.absolute(hs)) / np.sum(np.absolute(kx*hs)) #according to Duran(2007)
        
        # 1.2 Inner layer height
        l = 1.0
        for i in range(5):
            l_aux = 2. * kappa**2 * L / np.log(l/self.z0)
            l = l_aux
        
        # 1.3 Middle layer height
        hm = 1.0
        for i in range(5):
            hm_aux = L / np.sqrt(np.log(hm/self.z0))
            hm = hm_aux
        
        # 1.4 non-dimensional velocity
        ul = np.log(l/self.z0) / np.log(hm/self.z0)

        # 1.5 extra arryas in Fourier space
        k = np.sqrt(kx**2 + ky**2)
        sigma = np.sqrt(1j * L * kx * self.z0 / l)
        
        # 2. Shear stress perturbation
        #-----------------------------
        # According to Duran (2007)
        dtaux_t = hs * kx**2 / k * 2 / ul**2 * \
                  (-1 + (2 * np.log(l/self.z0) + k**2/kx**2) * sigma * \
                   scipy.special.kv(1, 2 * sigma) / scipy.special.kv(0, 2 * sigma))
        dtauy_t = hs * kx * ky / k * 2 / ul**2 * \
                  2 * np.sqrt(2) * sigma * scipy.special.kv(1, 2 * np.sqrt(2) * sigma)
        
        self.cgrid['dtaux'] = np.real(np.fft.ifft2(dtaux_t))
        self.cgrid['dtauy'] = np.real(np.fft.ifft2(dtauy_t))
        
    def set_computational_grid(self):
        '''Define computational grid
        
        The computational grid is square with dimensions equal to the
        diagonal of the bounding box of the input grid, plus twice the
        buffer width.

        '''

        gi = self.igrid
                
        # grid center
        x0, y0 = np.mean(gi['x']), np.mean(gi['y'])

        # grid size
        self.D = np.sqrt((gi['x'].max() - gi['x'].min())**2 +
                        (gi['y'].max() - gi['y'].min())**2) + 2 * self.buffer_width

        # determine equidistant, square grid
        xc, yc = self.get_exact_grid(x0 - self.D/2., x0 + self.D/2.,
                                    y0 - self.D/2., y0 + self.D/2.,
                                    self.cgrid['dx'], self.cgrid['dy'])

        self.x0 = x0
        self.y0 = y0
        self.cgrid['xi'] = xc
        self.cgrid['yi'] = yc
        
        # NO INTERPOLATION
        
        # self.x0 = x0
        # self.y0 = y0
        # self.cgrid['xi'] = g['x']
        # self.cgrid['yi'] = g['y']
        
        
    def get_sigmoid(self, x):
        '''Get sigmoid function value
        
        Get bed level multiplication factor in buffer area based on
        buffer specificationa and distance to input grid boundary.
        
        Parameters
        ----------
        x : float or numpy.ndarray
            Distance(s) to input grid boundary
        
        Returns
        -------
        float or numpy.ndarray
            Bed level multiplication factor (z = factor * z_boundary)

        '''
            
        return 1. / (1. + np.exp(-(self.buffer_width-x) / self.buffer_relaxation))
        

    def filter_highfrequenies(self, kx, ky, hs, nfilter=(1, 2), p=.01):
        '''Filter high frequencies from a 2D spectrum

        A logistic sigmoid filter is used to taper higher frequencies
        from the 2D spectrum. The range over which the sigmoid runs
        from 0 to 1 with a precision ``p`` is given by the 2-tuple
        ``nfilter``. The range is defined as wavenumbers in terms of
        gridcells, i.e. a value 1 corresponds to a wave with length
        ``dx``.

        Parameters
        ----------
        kx : numpy.ndarray
            Wavenumbers in x-direction
        ky : numpy.ndarray
            Wavenumbers in y-direction
        hs : numpy.ndarray
            2D spectrum
        nfilter : 2-tuple
            Wavenumber range used for logistic sigmoid filter
        p : float
            Precision of sigmoid range definition

        Returns
        -------
        hs : numpy.ndarray
            Filtered 2D spectrum

        '''
        
        if nfilter is not None:
            n1 = np.min(nfilter)
            n2 = np.max(nfilter)
            px = 2 * np.pi / self.cgrid['dx'] / np.abs(kx)
            py = 2 * np.pi / self.cgrid['dy'] / np.abs(ky)
            s1 =  n1 / np.log(1. / p - 1.)
            s2 = -n2 / np.log(1. / (1.- p) - 1.)
            f1 = 1. / (1. + np.exp(-(px + n1 - n2) / s1))
            f2 = 1. / (1. + np.exp(-(py + n1 - n2) / s2))
            hs *= f1 * f2

        return hs  
    
                                                                                                                                
    def plot(self, ax=None, cmap='Reds', stride=10, computational_grid=False, **kwargs):
        '''Plot wind shear perturbation
            
        Parameters
        ----------
        ax : matplotlib.pyplot.Axes, optional
            Axes to plot onto
        cmap : matplotlib.cm.Colormap or string, optional
            Colormap for topography (default: Reds)
        stride : int, optional
            Stride to apply to wind shear vectors (default: 10)
        computational_grid : bool, optional
            Plot on computational grid rather than input grid
            (default: False)
        kwargs : dict
            Additional arguments to :func:`matplotlib.pyplot.quiver`
            
        Returns
        -------
        ax : matplotlib.pyplot.Axes
            Axes used for plotting

        '''
        
        d = stride
        
        if ax is None:
            fig, ax = subplots()
        
        if computational_grid:
            g = self.cgrid
        else:
            g = self.igrid
        
        ax.pcolormesh(g['x'], g['y'], g['z'], cmap=cmap)
        ax.quiver(g['x'][::d,::d], g['y'][::d,::d], 
                  g['dtaux'][::d,::d], g['dtauy'][::d,::d], **kwargs)
                  
        if computational_grid:
            ax.plot(self.get_borders(self.igrid['x']),
                    self.get_borders(self.igrid['y']), '-k')
                  
        return ax


    @staticmethod
    def get_exact_grid(xmin, xmax, ymin, ymax, dx, dy):
        '''Returns a grid with given gridsizes approximately within given bounding box'''
        
        x = np.arange(np.floor(xmin / dx) * dx,
                      np.ceil(xmax / dx) * dx, dx)
        y = np.arange(np.floor(ymin / dy) * dy,
                      np.ceil(ymax / dy) * dy, dy)
        x, y = np.meshgrid(x, y)
                      
        return x, y
    
    
    @staticmethod
    def get_borders(x):
        '''Returns borders of a grid as one-dimensional array'''
        
        return np.concatenate((x[0,:].T, 
                               x[1:-1,-1], 
                               x[-1,::-1].T, 
                               x[-1:1:-1,0],
                               x[0,:1]), axis=0)
    
    
    @staticmethod
    def rotate(x, y, alpha, origin=(0,0)):
        '''Rotate a matrix over given angle around given origin'''
        
        xr = x - origin[0]
        yr = y - origin[1]
        
        a = alpha / 180. * np.pi
        
        R = np.asmatrix([[np.cos(a), -np.sin(a)],
                         [np.sin(a),  np.cos(a)]])
        
        xy = np.concatenate((xr.reshape((-1,1)), 
                             yr.reshape((-1,1))), axis=1) * R
                         
        return (np.asarray(xy[:,0].reshape(x.shape) + origin[0]),
                np.asarray(xy[:,1].reshape(y.shape) + origin[1]))
    
    
    def interpolate(self, x, y, z, xi, yi):
        '''Interpolate a grid onto another grid'''
        
        xy = np.concatenate((y.reshape((-1,1)),
                             x.reshape((-1,1))), axis=1)
        
        xyi = np.concatenate((yi.reshape((-1,1)),
                              xi.reshape((-1,1))), axis=1)  
        
        if self.istransect:
            zi = np.interp(xi.flatten(), x.flatten(), z.flatten()).reshape(xi.shape)
        else:
#            zi = scipy.interpolate.griddata(xy, z.reshape((-1,1)), xyi, method='cubic').reshape(xi.shape)
            inter = scipy.interpolate.RegularGridInterpolator((y[:,0], x[0,:]), z, bounds_error = False, fill_value = 0.)
            zi = inter(xyi).reshape(xi.shape)              
                     
        return zi


    @staticmethod
    def interpolate_projected_point(a, b, p):
        '''Project point to line segment and return distance and interpolated value
        
        Parameters
        ----------
        a : iterable
            Start vector for line segment
        b : iterable
            End vector for line segment
        p : iterable
            Point vector to be projected
            
        Returns
        -------
        d : numpy.ndarray
            Distance from point p to projected point q
        z : float
            Interpolated value at projected point q
            
        '''
        
        a = np.asarray(a)
        b = np.asarray(b)
        p = np.asarray(p)
        
        ab = b[:-1]-a[:-1]                     # line segment
        ab2 = np.dot(ab, ab)                   # length of line segment squared
        
        if ab2 > 0.:
            ap = p[:-1]-a[:-1]
            t = np.dot(ap, ab) / ab2           # location of projected point along line segment as fraction of its length
            if t >= 0. and t <= 1.:
                q = a[:-1] + t * ab            # projected point
                pq = p-q                       # vector from original to projected point
                d = np.sqrt(np.dot(pq, pq))    # distance from original to projected point
                z = a[-1] * (1.-t) + b[-1] * t # linearly interpolated height of projected point
                return d, z
            
        return
    
    def separation(self): 
        
        # Initialize grid and bed dimensions
        
        g = self.cgrid
         
        nx = len(g['z'][1])
        ny = len(g['z'][0])
        x = g['x']
        y = g['y']
        z = g['z']
        
        # Initialize arrays
    
        dzx  = np.zeros(g['z'].shape)
        dzy  = np.zeros(g['z'].shape)
        dz = np.zeros(g['z'].shape)
        stall = np.zeros(g['z'].shape)
        bubble = np.zeros(g['z'].shape)
        zsep = np.zeros(g['z'].shape) # total separation bubble
        zsep0 = np.zeros(g['z'].shape) # zero-order separation bubble
        zsep1 = np.zeros(g['z'].shape) # first-order separation bubble
        dzdx0 = np.zeros(g['z'].shape)
#        zsep2 = np.zeros(g['z'].shape) # separation bubble after cutting dune profile
        
#        zmin = np.zeros(ny)
        
        slope = 0.2 #np.deg2rad(11.) #p['M_dSlope'] #0.25
        
#        zfft = np.zeros(g['z'].shape, dtype=np.complex)
        k = np.array(range(0,nx))

        # Compute angle of separation
        dx = g['dx']
#
        dz[:,:-1] = np.rad2deg(np.arctan((z[:,1:]-z[:,:-1])/dx))
        dz[:,-1] = dz[:,-2]

        # dz = np.hypot(dzx,dzy)
        
        # Determine location of separation bubbles
        # stall += np.logical_and(dzx < 0, np.rad2deg(np.arctan(dz)) > 30.) #p['M_sep'] # 30

        stall += np.logical_and(abs(dz) > 20., dz < 0) # p['M_sep'] # 30
        
        stall[1:-1,:] += np.logical_and(stall[1:-1,:]==0, stall[:-2,:]>0, stall[2:,:]>0)
        stall[:,1:-1] += np.logical_and(stall[:,1:-1]==0, stall[:,:-2]>0, stall[:,2:]>0) 
        
        bubble[:,:-1] = np.logical_and(stall[:,:-1] == 0, stall[:,1:] > 0) # define bubble



        
        # Shift bubble n cells back
        n = 2
        bubble[:,:-n] = bubble[:,n:]
        bubble[:,:n] = 0
        
        bubble = bubble.astype(int)


        
        # Count separation bubbles
        
        n = np.sum(bubble)
        bubble_n = np.asarray(np.where(bubble == True)).T
        
        zfft = np.zeros((ny, nx), dtype=np.complex)
        
        for k in range(0, n):
            
            i = bubble_n[k,1]
            j = bubble_n[k,0]

            ix_neg = (dz[j, i:] >= 0)

            if np.sum(ix_neg) == 0:
                hb = z[j,i]
            else:
                hb = z[j,i] - z[j,i+np.where(ix_neg)[0][0]]

                # print('-----------------')
                # print(g['x'][j,i])
                # print(g['y'][j, i])
                # print(z[j,i])
                #
                # plt.pcolormesh(g['x'], g['y'], dz)
                # plt.colorbar()
                # plt.show()

            # Walk through all separation bubbles and determine polynoms
        
            dzdx0 = (z[j,i]-z[j,i-2])/(2*dx)
        
            a = dzdx0 / slope
            # l = np.minimum(np.maximum((1.5 * z[j,i] / slope) * (1 + a*0.25 + 0.125*a**2),.1),200.)
            l = np.minimum(np.maximum((1.5 * hb / slope) * (1 + a * 0.25 + 0.125 * a ** 2), .1), 200.)
            
            a2 = -3 * hb/l**2 - 2 * dzdx0 / l
            a3 =  2 * hb/l**3 +     dzdx0 / l**2
          
            i_max = min(i+int(l/dx),int(nx-1))
            # xs = x[j,i:i_max] - x[j,i]
            xs = np.arange(0, (i_max - i)*dx, dx)

            
            zsep0[j,i:i_max] = ((a3 * xs + a2) * xs + dzdx0) * xs + z[j,i]
            
            # First order Filter
#            
#            Cut = 1.5
#            dk = 2.0 * np.pi / (np.max(g['x']))
#            zfft[j,:] = np.fft.fft(zsep0[j,:])
#            zfft[j,:] *= np.exp(-(dk*k*g['dx'])**2./(2.*Cut**2.))
#            zsep0[j,:] = np.real(np.fft.ifft(zfft[j,:]))
#            
#            # First order
#            
#            dzdx1 = (zsep0[j,i+1]-zsep0[j,i])/(g['dx'])
#            
#            a = dzdx1 / slope
#            l = np.minimum(np.maximum((1.5 * z[j,i] / slope) * (1 + a*0.25 + 0.125*a**2),.1),200.)
#            
#            a2 = -3 * z[j,i]/l**2 - 2 * dzdx1 / l
#            a3 =  2 * z[j,i]/l**3 +     dzdx1 / l**2
#          
#            i_max = min(i+int(l/g['dx']),int(nx-1))
#            xs = x[j,i:i_max] - x[j,i]
#            
#            zsep1[j,i:i_max] = ((a3 * xs + a2) * xs + dzdx1) * xs + z[j,i]
            
            # Separation bubble cuts dune profile
            
#            cut_list = np.logical_and(zsep0[j, i:i_max] >= z[j,i:i_max], zsep0[j, i+1:i_max+1] < z[j,i:i_max])  
#            cut_list = cut_list.astype(int)
#            
#            if np.sum(cut_list[5:int(l[j,i]/g['dx'])]) > 0:
#            
#                i_cut_list = np.asarray(np.where(cut_list[5:] == True)).T
#                i_cut = int(i_cut_list[0]) + i
#    
#                dzdx1 = (z[j,i_cut] - z[j,i_cut-2])/(2.*g['dx'])
#                l_cut = x[j,i_cut] - x[j,i]
#                
#                a2_cut = (3.*z[j,i_cut]-dzdx1*l_cut-2.*dzdx0[j,i]*l_cut-3.*z[j,i])/(l_cut**2.)
#                a3_cut = (dzdx1*l_cut-2*z[j,i_cut]+dzdx0[j,i]*l_cut+2*z[j,i])/(l_cut**3.)
#                
#                i_max = min(i+int(l_cut/g['dx']),int(nx-1))
#                xs = x[j,i:i_max] - x[j,i]
#                zsep1[j,i:i_max] = ((a3_cut*xs + a2_cut) * xs + dzdx0[j,i])*xs + z[j,i]
#                zsep[j,i:i_max] = np.maximum(zsep[j,i:i_max], zsep1[j,i:i_max])
#            else:
            zsep[j,i:i_max] = np.maximum(zsep[j,i:i_max], zsep0[j,i:i_max]) # zsep1

        # plt.pcolormesh(g['x'], g['y'], zsep)
        # plt.colorbar()
        # plt.show()
        #
        # plt.plot(x[150,:], zsep[150,:])
        # plt.show()

        return zsep
    
    def separation_shear(self, dzsep):
        
        m_tau_sepbub = .05
        slope = 0.2 #np.tan(np.deg2rad(34.)) #Mcr_dyn
        delta = 1./(slope*m_tau_sepbub)
        
        zsepdelta = np.minimum(np.maximum(1. - delta * dzsep, 0.),1.)
        
        self.cgrid['taux'] *= zsepdelta
        self.cgrid['tauy'] *= zsepdelta
        
        return
    

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy.ma as ma
import warnings
import pyasdf
from pyproj import Geod, Proj
import obspy
from mpl_toolkits.basemap import Basemap
import colormaps
geodist = Geod(ellps='WGS84')
geoproj = Proj(proj='utm',zone=10,ellps='WGS84')

class vmodel(object):
    """
    An object to handle input model for SPECFEM2D.
    =======================================================================================================
    Parameters:
    xmin, xmax, zmin, zmax - bound of study region
    Nx, Nz                 - element number in x, z ([number of grid points]=[element number]*[lpd+1])
    Vp, Vs, Rho            - background Vp/Vs/density
    lpd                    - Lagrange polynomial degree
    plotflag               - whether to store numpy arrays for plotting purpose
    
    Notes:
    Nx, Nz is completely different from ni, nj in SW4. In SW4, ni, nj, nk are number of grid points,
    while Nx, Nz here are element numbers, which can be analogically regarded as block numbers.  
    =======================================================================================================
    """
    def __init__(self, xmin, xmax, Nx, zmin, zmax, Nz, Vp=4000., Vs=3000., Rho=2600., lpd=4, plotflag=True):
        self.xmin=xmin
        self.xmax=xmax
        self.Nx=Nx
        self.zmin=zmin
        self.zmax=zmax
        self.Nz=Nz
        self.Vp=Vp
        self.Vs=Vs
        self.Rho=Rho
        self.lpd=lpd
        self.dx=(xmax-xmin)/Nx
        self.dz=(zmax-zmin)/Nz
        self.XelementArr=self.xmin+np.arange(Nx)*self.dx+self.dx
        self.ZelementArr=self.zmin+np.arange(Nz)*self.dz+self.dz
        self.get_GLL()
        # Generate X grid points given element number and Lagrange polynomial degree
        A=np.ones([(lpd+1)*(lpd+1),Nx])
        repeatknots=np.tile(self.knots, (lpd+1))
        self.XArr=(A*self.XelementArr).T + (0.5+0.5*repeatknots*A.T)*self.dx - self.dx
        self.XArr=self.XArr.reshape( (lpd+1)*(lpd+1)*Nx )
        self.XArr=np.tile(self.XArr, Nz)
        # Generate X grid points given element number and Lagrange polynomial degree
        A=np.ones([lpd+1,Nz])
        B=(A*self.ZelementArr).T
        C=(0.5+0.5*A.T*self.knots)*self.dz - self.dz
        D=B+C
        E=np.tile(D, Nx)
        self.ZArr=np.repeat(E,lpd+1)
        # Generate background velocity and density parameters
        self.VpArr=np.ones( [( (lpd+1) * Nz )*( (lpd+1) * Nx )] )*Vp
        self.VsArr=np.ones( [( (lpd+1) * Nz )*( (lpd+1) * Nx )] )*Vs
        self.RhoArr=np.ones( [( (lpd+1) * Nz )*( (lpd+1) * Nx )] )*Rho
        self.plotflag=plotflag
        if plotflag:
        	A=np.ones([lpd+1, Nx])
        	A=self.dx*(0.5+0.5*(self.knots*A.T))
        	B=np.ones([lpd+1, Nx])
        	B=(B*self.XelementArr).T-self.dx
        	Xarr=A+B
        	
        	AA=np.ones([lpd+1, Nz])
        	AA=self.dz*(0.5+0.5*(self.knots*AA.T))
        	BB=np.ones([lpd+1, Nz])
        	BB=(BB*self.ZelementArr).T-self.dz
        	Zarr=AA+BB
        	
        	self.XArrPlot, self.ZArrPlot=np.meshgrid(Xarr,Zarr)
        	self.VpArrPlot=np.ones( [( (lpd+1) * Nz ), ( (lpd+1) * Nx )] )*Vp
        	self.VsArrPlot=np.ones( [( (lpd+1) * Nz ), ( (lpd+1) * Nx )] )*Vs
        	self.RhoArrPlot=np.ones( [( (lpd+1) * Nz ), ( (lpd+1) * Nx )] )*Rho
        self.regular=True
        return
    
    def read(self, infname):
        """
        Read txt velocity model file
        """
        self.regular=False
        InArr=np.loadtxt(infname)
        self.XArr=InArr[:,0]
        self.ZArr=InArr[:,1]
        self.RhoArr=InArr[:,2]
        self.VpArr=InArr[:,3]
        self.VsArr=InArr[:,4]
        self.Nin=self.RhoArr.size
        self.xmin=self.XArr.min()
        self.xmax=self.XArr.max()
        self.zmin=self.ZArr.min()
        self.zmax=self.ZArr.max()
        return
    
    def setbackground(self, vs=None, vp=None, rho=None):
        if vs !=None:
            self.VsArr[:]=vs
        if vp !=None:
            self.VpArr[:]=vp
        if rho !=None:
            self.RhoArr[:]=rho
        return
        
    
    
    def BlockHomoAnomaly(self, Xmin, Xmax, Zmin, Zmax, va=None, dv=None):
        """
        Inplement block anomaly in the model for Vs
        ============================================================================================
        Input Parameters:
        Xmin, Xmax, Zmin, Zmax - defines the bound
        va                     - anomalous velocity
        dv                     - velocity anomaly in percentage(default is None, which means use va)
        ============================================================================================
        """
        Xindex=(self.XArr>=Xmin)*(self.XArr<=Xmax)
        Zindex=(self.ZArr>=Zmin)*(self.ZArr<=Zmax)
        Index=Xindex*Zindex
        if dv!=None:
            self.VsArr[Index]=self.VsArr[Index]*(1+dv)
        else:
            self.VsArr[Index]=va
        if self.plotflag==True:
            Xindex=(self.XArrPlot>=Xmin)*(self.XArrPlot<=Xmax)
            Zindex=(self.ZArrPlot>=Zmin)*(self.ZArrPlot<=Zmax)
            Index=Xindex*Zindex
            if dv!=None:
                self.VsArrPlot[Index]=self.VsArrPlot[Index]*(1+dv)
            else:
                self.VsArrPlot[Index]=va
        return
    
    def CircleHomoAnomaly(self, Xc, Zc, R, va=None, dv=None):
        """
        Inplement circle anomaly in the model for Vs
        =============================================================================
        Input Parameters:
        Xc, Zc  - center of the circle
        R       - radius
        va      - anomalous velocity
        dv      - velocity anomaly in percentage(default is None, which means use va)
        =============================================================================
        """
        print 'Adding homo circle anomaly Xc=', Xc,' Zc=', Zc, ' R=',R
        dArr = np.sqrt( (self.XArr-Xc)**2 + (self.ZArr-Zc)**2)
        Index = dArr <= R
        if dv!=None:
            self.VsArr[Index]=self.VsArr[Index]*(1+dv)
        else:
            self.VsArr[Index]=va
        if self.plotflag==True:
            dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2)
            Index = dArr <= R
            if dv!=None:
                self.VsArrPlot[Index] = self.VsArrPlot[Index]*(1+dv)
            else:
                self.VsArrPlot[Index]=va
        return
    
    def CircleLinearAnomaly(self, Xc, Zc, R, va, dv=None):
        """
        Inplement circle anomaly with linear change towards center in the model for Vs
        Assuming the background Vs is homogeneous
        =============================================================================
        Input Parameters:
        Xc, Zc  - center of the circle
        R       - radius
        va      - anomalous velocity
        dv      - velocity anomaly in percentage(default is None, which means use va)
        =============================================================================
        """
        dArr = np.sqrt( (self.XArr-Xc)**2 + (self.ZArr-Zc)**2)
        if dv==None:
            dva = va - self.Vs
        else:
            dva = self.Vs*dv
        delD = R - dArr
        IndexIn = (delD >=0)
        # self.VsArr = 0.5 * (np.sign(delD) + 1) * delD/R * dva + self.VsArr
        self.VsArr = IndexIn * delD/R * dva + self.VsArr
        if self.plotflag==True:
            dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2)
            delD = R - dArr
            IndexIn = delD >=0
            # self.VsArrPlot = 0.5 * (np.sign(delD) + 1) * delD/R * dva + self.VsArrPlot
            self.VsArrPlot = IndexIn * delD/R * dva + self.VsArrPlot
        return
    
    def CircleCosineAnomaly(self, Xc, Zc, R, va=None, dv=None):
        """
        Inplement circle anomaly with cosine change towards center in the model for Vs
        Assuming the background Vs is homogeneous
        =============================================================================
        Input Parameters:
        Xc, Zc  - center of the circle
        R       - radius
        va      - anomalous velocity
        dv      - velocity anomaly in percentage(default is None, which means use va)
        =============================================================================
        """
        dArr = np.sqrt( (self.XArr-Xc)**2 + (self.ZArr-Zc)**2)
        if va !=None: dva = va - self.Vs
        else: dva = self.Vs*dv
        delD = R - dArr
        IndexIn = (delD >=0)
        # self.VsArr = 0.5 * (np.sign(delD) + 1) * ( 1+np.cos( np.pi* dArr / R ) ) * dva + self.VsArr
        self.VsArr = IndexIn * ( 1+np.cos( np.pi* dArr / R ) )/2. * dva + self.VsArr
        if self.plotflag==True:
            dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2)
            delD = R - dArr
            IndexIn = (delD >=0)
            # self.VsArrPlot = 0.5 * (np.sign(delD) + 1) * ( 1+np.cos( np.pi* dArr / R ) ) * dva + self.VsArrPlot
            self.VsArrPlot = IndexIn * ( 1+np.cos( np.pi* dArr / R ) )/2. * dva + self.VsArrPlot
        return
    
    def RingHomoAnomaly(self, Xc, Zc, Rmax, Rmin, va, dv=None):
        """
        Inplement ring anomaly in the model for Vs
        ================================================================================
        Input Parameters:
        Xc, Zc     - center of the circle
        Rmax, Rmin - radius max/min
        R          - radius
        va         - anomalous velocity
        dv         - velocity anomaly in percentage(default is None, which means use va)
        ================================================================================
        """
        if Rmin < self.dx:
            self.CircleHomoAnomaly(Xc=Xc, Zc=Zc, R=Rmax, va=va, dv=dv)
            return
        print 'Adding homo ring anomaly Xc=', Xc,' Zc=', Zc, ' Rmax=',Rmax, ' Rmin=', Rmin
        dArr = np.sqrt( (self.XArr-Xc)**2 + (self.ZArr-Zc)**2)
        Index = (dArr <= Rmax) * (dArr > Rmin) 
        if dv!=None: self.VsArr[Index]=self.VsArr[Index]*(1+dv)
        else: self.VsArr[Index]=va
        if self.plotflag==True:
            dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2)
            Index = (dArr <= Rmax) * (dArr > Rmin) 
            if dv!=None: self.VsArrPlot[Index] = self.VsArrPlot[Index]*(1+dv)
            else: self.VsArrPlot[Index]=va
        return
    
    def LensHomoModel(self, z0, x0, d1, d2, r1, r2, zmin=None, zmax=None, va=None, dv=None):
        if va==None and dv ==None:
            raise ValueError('va or dv need to be specified')
        x1      = x0-d1+r1
        x2      = x0+d2-r2
        L1      = np.sqrt(r1**2-(r1-d1)**2)
        L2      = np.sqrt(r2**2-(r2-d2)**2)
        zmax12  = z0+min(L1,L2)
        zmin12  = z0-min(L1,L2)
        if zmin!=None: zmin = max(zmin12, zmin)
        else: zmin  = zmin12
        if zmax!=None: zmax = min(zmax12, zmax)
        else: zmax  = zmax12
        zmin    = max(self.zmin, zmin)
        zmax    = min(self.zmax, zmax)
        print 'zmax =',zmax, 'zmin =', zmin
        ##################################################################
        Zindex  = (self.ZArr>=zmin)*(self.ZArr<=zmax)
        D1Arr   = np.sqrt( (self.XArr-x1)**2 + (self.ZArr-z0)**2)
        D2Arr   = np.sqrt( (self.XArr-x2)**2 + (self.ZArr-z0)**2)
        R1index = D1Arr<r1; R2index = D2Arr<r2
        Index   = R1index*R2index*Zindex
        if dv!=None: self.VsArr[Index]=self.VsArr[Index]*(1+dv)
        else: self.VsArr[Index]=va
        if self.plotflag==True:
            d1Arr = np.sqrt( (self.XArrPlot-x1)**2 + (self.ZArrPlot-z0)**2)
            d2Arr = np.sqrt( (self.XArrPlot-x2)**2 + (self.ZArrPlot-z0)**2)
            Index = (d1Arr < r1)*(d2Arr < r2)*(self.ZArrPlot>=zmin)*(self.ZArrPlot<=zmax)
            if dv!=None:
                self.VsArrPlot[Index] = self.VsArrPlot[Index]*(1+dv)
            else:
                self.VsArrPlot[Index]=va
        return
        
    def LensHomoModel_v2(self, z0, x0, d1, d2, r1, r2, zmin=None, zmax=None, va=None, dv=None):
        if va==None and dv ==None:
            raise ValueError('va or dv need to be specified')
        z1      = z0-d1+r1
        z2      = z0+d2-r2
        L1      = np.sqrt(r1**2-(r1-d1)**2)
        L2      = np.sqrt(r2**2-(r2-d2)**2)
        zmax12  = z0+min(L1,L2)
        zmin12  = z0-min(L1,L2)
        if zmin!=None: zmin = max(zmin12, zmin)
        else: zmin  = zmin12
        if zmax!=None: zmax = min(zmax12, zmax)
        else: zmax  = zmax12
        zmin    = max(self.zmin, zmin)
        zmax    = min(self.zmax, zmax)
        print 'zmax =',zmax, 'zmin =', zmin
        ##################################################################
        Zindex  = (self.ZArr>=zmin)*(self.ZArr<=zmax)
        D1Arr   = np.sqrt( (self.XArr-x0)**2 + (self.ZArr-z1)**2)
        D2Arr   = np.sqrt( (self.XArr-x0)**2 + (self.ZArr-z2)**2)
        R1index = D1Arr<r1; R2index = D2Arr<r2
        Index   = R1index*R2index*Zindex
        if dv!=None: self.VsArr[Index]=self.VsArr[Index]*(1+dv)
        else: self.VsArr[Index]=va
        if self.plotflag==True:
            d1Arr = np.sqrt( (self.XArrPlot-x0)**2 + (self.ZArrPlot-z1)**2)
            d2Arr = np.sqrt( (self.XArrPlot-x0)**2 + (self.ZArrPlot-z2)**2)
            Index = (d1Arr < r1)*(d2Arr < r2)*(self.ZArrPlot>=zmin)*(self.ZArrPlot<=zmax)
            if dv!=None:
                self.VsArrPlot[Index] = self.VsArrPlot[Index]*(1+dv)
            else:
                self.VsArrPlot[Index]=va
        return
    
    def HomoEllipse(self, Xc, Zc, a, b, va=None, dv=None):
        if va==None and dv ==None:
            raise ValueError('va or dv need to be specified')
        print 'Adding homo ellipse anomaly Xc =', Xc,' Zc =', Zc, ' a =',a, 'b =',b
        unitArr = np.sqrt( (self.XArr-Xc)**2/a**2 + (self.ZArr-Zc)**2/b**2)
        Index = unitArr <= 1
        if dv!=None:
            self.VsArr[Index]=self.VsArr[Index]*(1+dv)
        else:
            self.VsArr[Index]=va
        if self.plotflag==True:
            dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2)
            Index = dArr <= R
            if dv!=None:
                self.VsArrPlot[Index] = self.VsArrPlot[Index]*(1+dv)
            else:
                self.VsArrPlot[Index]=va
        
        
        z1      = z0-d1+r1
        z2      = z0+d2-r2
        L1      = np.sqrt(r1**2-(r1-d1)**2)
        L2      = np.sqrt(r2**2-(r2-d2)**2)
        zmax12  = z0+min(L1,L2)
        zmin12  = z0-min(L1,L2)
        if zmin!=None: zmin = max(zmin12, zmin)
        else: zmin  = zmin12
        if zmax!=None: zmax = min(zmax12, zmax)
        else: zmax  = zmax12
        zmin    = max(self.zmin, zmin)
        zmax    = min(self.zmax, zmax)
        print 'zmax =',zmax, 'zmin =', zmin
        ##################################################################
        Zindex  = (self.ZArr>=zmin)*(self.ZArr<=zmax)
        D1Arr   = np.sqrt( (self.XArr-x0)**2 + (self.ZArr-z1)**2)
        D2Arr   = np.sqrt( (self.XArr-x0)**2 + (self.ZArr-z2)**2)
        R1index = D1Arr<r1; R2index = D2Arr<r2
        Index   = R1index*R2index*Zindex
        if dv!=None: self.VsArr[Index]=self.VsArr[Index]*(1+dv)
        else: self.VsArr[Index]=va
        if self.plotflag==True:
            d1Arr = np.sqrt( (self.XArrPlot-x0)**2 + (self.ZArrPlot-z1)**2)
            d2Arr = np.sqrt( (self.XArrPlot-x0)**2 + (self.ZArrPlot-z2)**2)
            Index = (d1Arr < r1)*(d2Arr < r2)*(self.ZArrPlot>=zmin)*(self.ZArrPlot<=zmax)
            if dv!=None:
                self.VsArrPlot[Index] = self.VsArrPlot[Index]*(1+dv)
            else:
                self.VsArrPlot[Index]=va
        return
        
        
    
    def readPhv(self, x0, z0, infname, evlo=None, evla=None, lon0=None, lat0=None, dx=0.5):
        inArr   = np.loadtxt(infname)
        lons    = inArr[:,0]; lats=inArr[:,1]; phvArr=inArr[:,2]*1000.
        lon0    = lons.min(); lat0=lats.min(); lon1=lons.max(); lat1=lats.max()
        del_lons= lons-lon0; del_lats=lats-lat0; maxdlon=del_lons.max(); maxdlat=del_lats.max()
        # maxdistNS=obspy.geodetics.gps2dist_azimuth(lat0, lon0, lat0+maxdlat, lon0) # distance is in m
        # maxdistEW=obspy.geodetics.gps2dist_azimuth(lat0, lon0, lat0, lon0+maxdlon) # distance is in m
        # maxdist_y=maxdistNS/1000.; maxdist_x=maxdistEW/1000.
        m       = Basemap(projection='cea',llcrnrlat=lat0,urcrnrlat=lat1, llcrnrlon=lon0,urcrnrlon=lon1,resolution='c')
        xins, zins = m(lons, lats)
        try:
            evx, evz=m(evlo, evla)
            print 'Source location: ', evx, evz
        except:
            pass
        # if xins.max()>self.xmax-x0 or zins.max()>self.zmax-z0:
        #     raise ValueError('Input phase velocity map is too large!')
        Numb    = int(lons.size)
        radius, az, baz=obspy.geodetics.gps2dist_azimuth(0, 0, dx, 0)
        for i in xrange(Numb):
            # print i
            x   = xins[i]; z=zins[i]; lon=lons[i]; lat=lats[i]; va=phvArr[i]
            # xc=x0+x; zc=z0+z
            # self.CircleHomoAnomaly(Xc=xc, Zc=zc, R=radius, va=va)
            xmin= x0+x-radius; xmax=x0+x+radius; zmin=z0+z-radius; zmax=z0+z+radius
            self.BlockHomoAnomaly(Xmin=xmin, Xmax=xmax, Zmin=zmin, Zmax=zmax, va=va)
        return
    
    def ASDFmodel(self, infname, per=10., phgr=1, verbose=True, Dx=0., Dz=0.):
        """
        Read ASDF model
        =============================================================================
        Input Parameters:
        infname   - input file name
        per       - period
        phgr      - use phase(1) or group(2) velocity
        =============================================================================
        """
        dbase = pyasdf.ASDFDataSet(infname)
        if verbose==True:
            print dbase.auxiliary_data.Disp
        perArr = dbase.auxiliary_data.Disp.VP000.data.value[0, :]
        vArr=dbase.auxiliary_data.Disp.VP000.data.value[phgr, :]
        if not np.any(perArr==per):
            raise ValueError('Period ', per,' sec not in the theoretical dispersion curve !' )
        Vs0 = vArr[perArr==per]* 1000.
        self.setbackground(vs=Vs0) 
        if self.plotflag == True: self.VsArrPlot[:]=Vs0 
        for auxid in dbase.auxiliary_data.Disp.list()[1:]:
            perArr = dbase.auxiliary_data.Disp[auxid].data.value[0, :]
            vArr=dbase.auxiliary_data.Disp[auxid].data.value[phgr, :]
            Rmax=dbase.auxiliary_data.Disp[auxid].parameters['Rmax']
            Rmin=dbase.auxiliary_data.Disp[auxid].parameters['Rmin']
            x=dbase.auxiliary_data.Disp[auxid].parameters['x']
            y=dbase.auxiliary_data.Disp[auxid].parameters['y']
            Vs = vArr[perArr==per] * 1000.
            self.RingHomoAnomaly(Xc=x+Dx, Zc=y+Dz, Rmax=Rmax, Rmin=Rmin, va=Vs)
        return
            
    
    def write(self, outfname, dt=None, fc=None, freqfactor=2.5, C=0.35):
        """
        Write vmodel to txt file that can be read by SPECFEM2D
        Will check stability criteria if dt, fc are specified.
        =============================================================================
        Input Parameters:
        outfname   - output file name
        dt         - time step 
        fc         - center frequency
        freqfactor - fmin = fc * freqfactor
        C          - Courant number
        =============================================================================
        """
        if dt ==None and fc == None:
            warnings.warn('Skip input checker, may cause problem in simulation!', UserWarning, stacklevel=1)
        else:
            vmin, vmax = self.GetMinMaxV()
            checker=InputChecker(dt=dt, dx=self.dx, dz=self.dz, fc=fc, lpd=self.lpd, vmin=vmin, vmax=vmax)
            checker.Check(freqfactor=freqfactor, C=C)
        N=self.XArr.size
        OutArr=np.append(self.XArr, self.ZArr)
        OutArr=np.append(OutArr, self.RhoArr)
        OutArr=np.append(OutArr, self.VpArr)
        OutArr=np.append(OutArr, self.VsArr)
        OutArr=OutArr.reshape( 5, N )
        OutArr=OutArr.T
        np.savetxt(outfname, OutArr, fmt='  %.6e')
        return
    
    def get_GLL(self):
        """
        Set Gauss-Lobatto-Legendre(GLL) points for a given Lagrange polynomial degree.
        To construct a polynomial of degree n passing through n+1 data points. 
        """
        if self.lpd == 2:
            knots = np.array([-1.0, 0.0, 1.0])
        elif self.lpd == 3:
            knots = np.array([-1.0, -0.4472135954999579, 0.4472135954999579, 1.0])
        elif self.lpd == 4:
            knots = np.array([-1.0, -0.6546536707079772, 0.0, 0.6546536707079772, 1.0])
            # knots = np.array([-1.0, -1./7.*np.sqrt(21.), 0.0, 1./7.*np.sqrt(21.), 1.0])
        elif self.lpd == 5:
            knots = np.array([-1.0, -0.7650553239294647, -0.2852315164806451, 0.2852315164806451, 0.7650553239294647, 1.0])
        elif self.lpd == 6:
            knots = np.array([-1.0, -0.8302238962785670, -0.4688487934707142, 0.0, 0.4688487934707142, 0.8302238962785670, 1.0])
        elif self.lpd == 7:
            knots = np.array([-1.0, -0.8717401485096066, -0.5917001814331423,\
                -0.2092992179024789, 0.2092992179024789, 0.5917001814331423, 0.8717401485096066, 1.0])
        self.knots=knots
        return
    
    def plot(self, ds=1000, unit='km', cmap='seismic_r', vmin=None, vmax=None, zsize=10):
        """Plot velocity model
        =============================================================================
        Input Parameters:
        ds              - grid spacing
        unit            - unit
        vmin, vmax      - vmin,vmax for colorbar
        =============================================================================
        """
        
        # XLength=self.xmax-self.xmin
        # ZLength=self.zmax-self.zmin
        # xsize=zsize*(XLength/ZLength)
        # fig = plt.figure(figsize=(xsize, zsize))
        if cmap=='ses3d':
            cmap = colormaps.make_colormap({0.0:[0.1,0.0,0.0], 0.2:[0.8,0.0,0.0], 0.3:[1.0,0.7,0.0],0.48:[0.92,0.92,0.92],
                0.5:[0.92,0.92,0.92], 0.52:[0.92,0.92,0.92], 0.7:[0.0,0.6,0.7], 0.8:[0.0,0.0,0.8], 1.0:[0.0,0.0,0.1]})
        if self.plotflag==False:
            raise ValueError('No plot array!')
        # plt.figure(figsize=(16,13))
        if self.regular==True:
            im=plt.pcolormesh(self.XArrPlot/ds, self.ZArrPlot/ds, self.VsArrPlot/ds, cmap=cmap, vmin=vmin, vmax=vmax)
        else:
            xi = np.linspace(self.xmin, self.xmax, self.Nx*10)
            zi = np.linspace(self.zmin, self.zmax, self.Nz*10)
            self.xi, self.zi = np.meshgrid(xi, zi)
            #-- Interpolating at the points in xi, yi
            self.vi = griddata(self.XArr, self.ZArr, self.VsArr, self.xi, self.zi, 'linear')
            im=plt.pcolormesh(self.xi/ds, self.zi/ds, ma.getdata(self.vi)/ds, cmap=cmap, vmin=vmin, vmax=vmax, shading='gouraud')
        ##########################################
        plt.plot(500., 1000 , 'r*', markersize=30, lw=3)
        # plt.plot( [0., 4000.], [1000, 1000] , 'b--', lw=3)
        # plt.plot( [500., 500.], [700., 1300.] , 'g-', lw=3)
        # plt.plot( [500., 3500.], [700., 700.] , 'g-', lw=3)
        # plt.plot( [500., 3500.], [1300., 1300.] , 'g-', lw=3)
        # plt.plot( [3500., 3500.], [700., 1300.] , 'g-', lw=3)
        # 
        # plt.plot( [0., 0.], [0., 2000.] , 'k-', lw=3)
        # plt.plot( [0., 4000.], [0., 0.] , 'k-', lw=3)
        # plt.plot( [4000., 4000.], [0., 2000.] , 'k-', lw=3)
        # plt.plot( [0., 4000.], [2000., 2000.] , 'k-', lw=3)
        ##########################################
        plt.xlabel('x('+unit+')', fontsize=35)
        plt.ylabel('z('+unit+')', fontsize=35)
        # plt.axis([self.xmin/ds, self.xmax/ds, self.zmin/ds, self.zmax/ds])
        plt.axis('scaled')
        cb=plt.colorbar(shrink=0.8)#, size="3%", pad='2%')
        cb.set_label('Vs (km/s)', fontsize=20, rotation=90)
        plt.yticks(fontsize=30)
        plt.xticks(fontsize=30)
        ########################
        # plt.ylim([-100, 2100])
        # plt.xlim([-100, 4100])
        ########################
        plt.show()
        return
    

    def GetMinMaxV(self):
        """
        Get minimum/maximum vs 
        """
        vmin=self.VsArr.min()
        vmax=self.VsArr.max()
        return vmin, vmax
    
class InputChecker(object):
    """
    An object to check stability condition given input parameters.
    =============================================================================
    Parameters:
    dt           - time step
    dx, dz       - element spacing 
    fc           - central frequency
    lpd          - Lagrange polynomial degree
    vmin, vmax   - minimum/maximum velocity
    =============================================================================
    """
    def __init__(self, dt, dx, dz, fc, lpd, vmin, vmax):
        self.dt=dt
        self.dx=dx
        self.dz=dz
        self.fc=fc
        self.lpd=lpd
        self.get_GLL()
        self.vmin=vmin
        self.vmax=vmax
        return
    
    def get_GLL(self):
        """
        Set Gauss-Lobatto-Legendre(GLL) points for a given Lagrange polynomial degree.
        To construct a polynomial of degree n passing through n+1 data points. 
        """
        if self.lpd == 2:
            knots = np.array([-1.0, 0.0, 1.0])
        elif self.lpd == 3:
            knots = np.array([-1.0, -0.4472135954999579, 0.4472135954999579, 1.0])
        elif self.lpd == 4:
            knots = np.array([-1.0, -0.6546536707079772, 0.0, 0.6546536707079772, 1.0])
        elif self.lpd == 5:
            knots = np.array([-1.0, -0.7650553239294647, -0.2852315164806451, 0.2852315164806451, 0.7650553239294647, 1.0])
        elif self.lpd == 6:
            knots = np.array([-1.0, -0.8302238962785670, -0.4688487934707142, 0.0, 0.4688487934707142, 0.8302238962785670, 1.0])
        elif self.lpd == 7:
            knots = np.array([-1.0, -0.8717401485096066, -0.5917001814331423,\
                -0.2092992179024789, 0.2092992179024789, 0.5917001814331423, 0.8717401485096066, 1.0])
        self.knots=knots
        return
    
    def CheckMinLambda(self, freqfactor=2.5):
        """
        Check grid spacing with wavelength minimum wavelength.
        ==============================================
        Input Parameters:
        freqfactor       - fmax = freqfactor*fc
        ==============================================
        """
        lambdamin=self.vmin/self.fc/freqfactor
        dxArr=self.dx*np.diff( (0.5+0.5*(self.knots)) )
        dzArr=self.dz*np.diff( (0.5+0.5*(self.knots)) )
        dsmax=max(dxArr.max(), dzArr.max())
        # dsmax=max(self.dx, self.dz)
        # Need checking! (in manual: threshold value is around 4.5 points
        # per wavelength in elastic media and 5.5 in acoustic media), 4.5 grid points OR 4.5 element points
        if dsmax * 4.5 > lambdamin:
            raise ValueError('Grid spacing is too large: '+ str(dsmax)+' for '+str(lambdamin)+ ' m')
        else:
            print 'Grid spacing:', str(dsmax),'m for',lambdamin, 'm'
        return
    
    def CheckCFLCondition(self, C=0.35):
        """
        Check Courant-Frieddrichs-Lewy stability condition
        ====================================================
        Input Parameters:
        C - Courant number (default = 0.35, normally 0.3~0.4)
        ====================================================
        """
        dxArr=self.dx*np.diff( (0.5+0.5*(self.knots)) )
        dzArr=self.dz*np.diff( (0.5+0.5*(self.knots)) )
        dsmin=min(dxArr.min(), dzArr.min())
        dtCFL=C*dsmin/self.vmax
        if self.dt > dtCFL:
            raise ValueError('Time step violates Courant-Frieddrichs-Lewy Condition: ', dt, dtCFL)
        else:
            print 'Time Step: ',self.dt,' s Required Time Step: ', dtCFL, 's'
        return
    
    def Check(self, freqfactor=2.5, C=0.35):
        """
        Check minimum wavelenght and Courant conditions
        """
        print '=========== Checking stability conditions ==========='
        self.CheckMinLambda(freqfactor=freqfactor)
        self.CheckCFLCondition(C=C)
        print '===========  Stability conditions checked  ==========='
        return



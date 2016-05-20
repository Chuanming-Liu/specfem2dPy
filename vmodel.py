import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy.ma as ma
    
class vmodel(object):
    """
    An object to handle input model for SPECFEM2D.
    --------------------------------------------------------------------------------------------------------------
    Parameters:
    xmin, xmax, zmin, zmax    - bound of study region
    Nx, Nz                                - element number in x, z
                                                ( [number of grid points] = [element number] * [lpd+1] )
    Vp, Vs, Rho                        - background Vp/Vs/density
    lpd                                     - Lagrange polynomial degree
    plotflag                              - whether to store numpy arrays for plotting purpose 
    --------------------------------------------------------------------------------------------------------------
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
        if plotflag==True:
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
    
    def BlockHomoAnomaly(self, Xmin, Xmax, Zmin, Zmax, va, dv=None):
        """
        Inplement block anomaly in the model for Vs
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        Xmin, Xmax, Zmin, Zmax - defines the bound
        va                                     - anomalous velocity
        dv                                     - velocity anomaly in percentage( default is None, which means use va )
        -----------------------------------------------------------------------------------------------------
        """
        Xlogic=(self.XArr>=Xmin)*(self.XArr<=Xmax)
        Zlogic=(self.ZArr>=Zmin)*(self.ZArr<=Zmax)
        Logic=Xlogic*Zlogic
        if dv!=None:
            self.VsArr[Logic]=self.VsArr[Logic]*(1+dv)
        else:
            self.VsArr[Logic]=va
        if self.plotflag==True:
            Xlogic=(self.XArrPlot>=Xmin)*(self.XArrPlot<=Xmax)
            Zlogic=(self.ZArrPlot>=Zmin)*(self.ZArrPlot<=Zmax)
            Logic=Xlogic*Zlogic
            if dv!=None:
                self.VsArrPlot[Logic]=self.VsArrPlot[Logic]*(1+dv)
            else:
                self.VsArrPlot[Logic]=va
        return
    
    def CircleHomoAnomaly(self, Xc, Zc, R, va, dv=None):
        """
        Inplement circle anomaly in the model for Vs
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        Xc, Zc      - the center of the circle
        R             - radius
        va           - anomalous velocity
        dv           - velocity anomaly in percentage( default is None, which means use va )
        -----------------------------------------------------------------------------------------------------
        """
        dArr = np.sqrt( (self.XArr-Xc)**2 + (self.ZArr-Zc)**2)
        Logic = dArr < R
        if dv!=None:
            self.VsArr[Logic]=self.VsArr[Logic]*(1+dv)
        else:
            self.VsArr[Logic]=va
        if self.plotflag==True:
            dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2)
            Logic = dArr < R
            if dv!=None:
                self.VsArrPlot[Logic] = self.VsArrPlot[Logic]*(1+dv)
            else:
                self.VsArrPlot[Logic]=va
        return
    
    def CircleLinearAnomaly(self, Xc, Zc, R, va, dv=None):
        """
        Inplement circle anomaly with linear change towards center in the model for Vs
        Assuming the background Vs is homogeneous
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        Xc, Zc      - the center of the circle
        R             - radius
        va           - anomalous velocity
        dv           - velocity anomaly in percentage( default is None, which means use va )
        -----------------------------------------------------------------------------------------------------
        """
        dArr = np.sqrt( (self.XArr-Xc)**2 + (self.ZArr-Zc)**2)
        if dv!=None:
            dva = va - self.Vs
        else:
            dva = self.Vs*dv
        delD = R - dArr 
        self.VsArr = 0.5 * (np.sign(delD) + 1) * delD/R * dva + self.VsArr
        if self.plotflag==True:
            dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2)
            delD = R - dArr
            self.VsArrPlot = 0.5 * (np.sign(delD) + 1) * delD/R * dva + self.VsArrPlot
        return
    
    def CircleCosineAnomaly(self, Xc, Zc, R, va, dv=None):
        """
        Inplement circle anomaly with cosine change towards center in the model for Vs
        Assuming the background Vs is homogeneous
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        Xc, Zc      - the center of the circle
        R             - radius
        va           - anomalous velocity
        dv           - velocity anomaly in percentage( default is None, which means use va )
        -----------------------------------------------------------------------------------------------------
        """
        dArr = np.sqrt( (self.XArr-Xc)**2 + (self.ZArr-Zc)**2)
        if dv!=None:
            dva = va - self.Vs
        else:
            dva = self.Vs*dv
        delD = R - dArr 
        self.VsArr = 0.5 * (np.sign(delD) + 1) * ( 1+np.cos( np.pi* dArr / R ) ) * dva + self.VsArr
        if self.plotflag==True:
            dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2)
            delD = R - dArr
            self.VsArrPlot = 0.5 * (np.sign(delD) + 1) * ( 1+np.cos( np.pi* dArr / R ) ) * dva + self.VsArrPlot
        return
    
    def write(self, outfname):
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
    
    def plot(self, ds=1000, unit='km', vmin=2.5, vmax=3.5):
        if self.plotflag==False:
            raise ValueError('No plot array!')
        plt.figure(figsize=(16,12))
        if self.regular==True:
            plt.pcolormesh(self.XArrPlot/ds, self.ZArrPlot/ds, self.VsArrPlot/ds, cmap='seismic_r', vmin=vmin, vmax=vmax)
        else:
            xi = np.linspace(self.xmin, self.xmax, self.Nx*10)
            zi = np.linspace(self.zmin, self.zmax, self.Nz*10)
            self.xi, self.zi = np.meshgrid(xi, zi)
            #-- Interpolating at the points in xi, yi
            self.vi = griddata(self.XArr, self.ZArr, self.VsArr, self.xi, self.zi, 'linear')
            plt.pcolormesh(self.xi/ds, self.zi/ds, ma.getdata(self.vi)/ds, cmap='seismic_r', vmin=vmin, vmax=vmax)
        # plt.plot( 320, 320 , 'y*', markersize=30)
        plt.xlabel('x('+unit+')', fontsize=30)
        plt.ylabel('z('+unit+')', fontsize=30)
        plt.colorbar()
        plt.axis([self.xmin/ds, self.xmax/ds, self.zmin/ds, self.zmax/ds])
        # plt.axis('scaled')
        plt.yticks(fontsize=20)
        plt.xticks(fontsize=20)
        plt.show()
        return
    
    def GetMinMaxV(self):
        vmin=self.VsArr.min()
        vmax=self.VsArr.max()
        return vmin, vmax
    
class InputChecker(object):
    """
    An object to check stability condition given input parameters.
    -----------------------------------------------------------------------------------------------------
    Parameters:
    dt                  - time step
    dx, dz            - grid spacing 
    fc                   - central frequency
    lpd                 - Lagrange polynomial degree
    vmin, vmax   - minimum/maximum velocity
    -----------------------------------------------------------------------------------------------------
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
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        freqfactor - fmax = freqfactor*fc
        -----------------------------------------------------------------------------------------------------
        """
        lambdamin=self.vmin/self.fc/freqfactor
        dxArr=self.dx*np.diff( (0.5+0.5*(self.knots)) )
        dzArr=self.dz*np.diff( (0.5+0.5*(self.knots)) )
        dsmax=max(dxArr.max(), dzArr.max())
        if dsmax * 10. > lambdamin:
            raise ValueError('Grid spacing is too large: ', str(dsmax),' for ',lambdamin, ' km')
        else:
            print 'Grid spacing: ', str(dsmax),'km for ',lambdamin, ' km'
        # if self.dx * 2. > lambdamin or self.dy * 2. > lambdamin:
        #     raise ValueError('Grid spacing is too large: ', str(dsmax),' for ',lambdamin, ' km')
        # else:
        #     print 'Grid spacing: ', str(dsmax),'km for ',lambdamin, ' km'
        return
    
    def CheckCFLCondition(self, C=0.35):
        """
        Check Courant-Frieddrichs-Lewy stability condition
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        C - Courant number (default = 0.35, normally 0.3~0.4)
        -----------------------------------------------------------------------------------------------------
        """
        dxArr=self.dx*np.diff( (0.5+0.5*(self.knots)) )
        dzArr=self.dz*np.diff( (0.5+0.5*(self.knots)) )
        dsmin=min(dxArr.min(), dzArr.min())
        dtCFL=C*dsmin/self.vmax
        if self.dt > dtCFL:
            raise ValueError('Time step violates Courant-Frieddrichs-Lewy Condition: ', dt, dtCFL)
        else:
            print 'Time Step: ',self.dt,' s Required Time Step: ', dtCFL, 's '
        return
    
    def Check(self, freqfactor=2.5, C=0.35):
        """
        Check minimum wavelenght and Courant conditions
        """
        print '-------------------- Checking stability conditions --------------------'
        self.CheckMinLambda(freqfactor=freqfactor)
        self.CheckCFLCondition(C=C)
        print '-------------------- Stability conditions checked  ---------------------'
        return



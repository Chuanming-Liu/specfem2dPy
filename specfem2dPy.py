import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy.ma as ma
import matplotlib.animation as animation

class StaInfo(object):
    """
    An object contains a station information several methods for station related analysis.
    -----------------------------------------------------------------------------------------------------
    General Parameters:
    stacode     - station name
    network     - network
    chan        - channels for analysis
    x,z     - position for station
    -----------------------------------------------------------------------------------------------------
    """
    def __init__(self, stacode=None, network='ME2D', x=None, z=None):

        self.stacode=stacode;
        self.network=network;
        self.x=x;
        self.z=z;
        return;
    
class StaLst(object):
    """
    An object contains a station list(a list of StaInfo object) information several methods for station list related analysis.
        stations: list of StaInfo
    """
    def __init__(self,stations=None):
        self.stations=[]
        if isinstance(stations, StaInfo):
            stations = [stations]
        if stations:
            self.stations.extend(stations)

    def __add__(self, other):
        """
        Add two StaLst with self += other.
        """
        if isinstance(other, StaInfo):
            other = StaLst([other])
        if not isinstance(other, StaLst):
            raise TypeError
        stations = self.stations + other.stations
        return self.__class__(stations=stations)

    def __len__(self):
        """
        Return the number of Traces in the StaLst object.
        """
        return len(self.stations)

    def __getitem__(self, index):
        """
        __getitem__ method of obspy.Stream objects.
        :return: Trace objects
        """
        if isinstance(index, slice):
            return self.__class__(stations=self.stations.__getitem__(index))
        else:
            return self.stations.__getitem__(index)

    def append(self, station):
        """
        Append a single StaInfo object to the current StaLst object.
        """
        if isinstance(station, StaInfo):
            self.stations.append(station)
        else:
            msg = 'Append only supports a single StaInfo object as an argument.'
            raise TypeError(msg)
        return self

    def ReadStaList(self, stafile):
        """
        Read Sation List from a txt file
        stacode network x z
        """
        f = open(stafile, 'r')
        Sta=[]
        for lines in f.readlines():
            lines=lines.split()
            stacode=lines[0]
            network=lines[1]
            x=float(lines[2])
            z=float(lines[3])
            # if len(lines)==5:
            #     try:
            #         ccflag=int(lines[3])
            #         network=lines[4]
            #     except ValueError:
            #         ccflag=int(lines[4])
            #         network=lines[3]
            # if len(lines)==4:
            #     try:
            #         ccflag=int(lines[3])
            #     except ValueError:
            #         network=lines[3]
            netsta=network+'.'+stacode
            if Sta.__contains__(netsta):
                index=Sta.index(netsta)
                if abs(self[index].lon-lon) >0.01 and abs(self[index].lat-lat) >0.01:
                    raise ValueError('Incompatible Station Location:' + netsta+' in Station List!')
                else:
                    print 'Warning: Repeated Station:' +netsta+' in Station List!'
                    continue
            Sta.append(netsta)
            self.append(StaInfo (stacode=stacode, network=network, x=x, z=z ))
            f.close()
        return
    
    def WriteStaList(self, stafile):
        """
        Read Sation List from a txt file
        stacode network x z
        """
        f = open(stafile, 'w')
        for sta in self.stations:
                tempstr='%s %s %1.5e %1.5e 0.0 0.0  \n' %( sta.stacode, sta.network, sta.x, sta.z )
                f.writelines(tempstr)
        f.close()
        return
    
    def AddSingle(self, x, z, stacode=None, network='MEM2D'):
        if stacode==None:
            stacode=str(int(x))+'S'+str(int(z));
        self.append(StaInfo (stacode=stacode, network=network, x=x, z=z ));
        return;
    
    def HomoStaLst(self, xmin, xmax, dx, zmin, zmax, dz, network='MEM2D'):
        Nx = int((xmax-xmin)/dx)+1;
        Nz = int((zmax-zmin)/dz)+1;
        for i in np.arange(Nx):
            for j in np.arange(Nz):
                x=xmin+dx*i;
                z=zmin+dz*j;
                stacode=str(int(i))+'S'+str(int(j));
                self.append(StaInfo (stacode=stacode, network=network, x=x, z=z ));
        return;
    
class VelocityModel(object):
    
    def __init__(self, xmin, xmax, Nx, zmin, zmax, Nz, Vp=4000., Vs=3500., Rho=2600., lpd=4, plotflag=True):
        self.xmin=xmin;
        self.xmax=xmax;
        self.Nx=Nx;
        self.zmin=zmin;
        self.zmax=zmax;
        self.Nz=Nz;
        self.Vp=Vp;
        self.Vs=Vs;
        self.Rho=Rho;
        self.lpd=lpd;
        self.dx=(xmax-xmin)/Nx;
        self.dz=(zmax-zmin)/Nz;
        self.XelementArr=self.xmin+np.arange(Nx)*self.dx+self.dx;
        self.ZelementArr=self.zmin+np.arange(Nz)*self.dz+self.dz;
        self.get_GLL();
        
        A=np.ones([(lpd+1)*(lpd+1),Nx]);
        repeatknots=np.tile(self.knots, (lpd+1));
        self.XArr=(A*self.XelementArr).T + (0.5+0.5*repeatknots*A.T)*self.dx - self.dx;
        self.XArr=self.XArr.reshape( (lpd+1)*(lpd+1)*Nx );
        self.XArr=np.tile(self.XArr, Nx)
    
        A=np.ones([lpd+1,Nx])
        B=(A*self.ZelementArr).T
        C=(0.5+0.5*A.T*self.knots)*self.dx - self.dx;
        D=B+C
        E=np.tile(D, Nx)
        self.ZArr=np.repeat(E,lpd+1)
        self.VpArr=np.ones( [( (lpd+1) * Nz )*( (lpd+1) * Nx )] )*Vp;
        self.VsArr=np.ones( [( (lpd+1) * Nz )*( (lpd+1) * Nx )] )*Vs;
        self.RhoArr=np.ones( [( (lpd+1) * Nz )*( (lpd+1) * Nx )] )*Rho;

        self.plotflag=plotflag;
        if plotflag==True:
        	A=np.ones([lpd+1, Nx]);
        	A=self.dx*(0.5+0.5*(self.knots*A.T));
        	B=np.ones([lpd+1, Nx]);
        	B=(B*self.XelementArr).T-self.dx;
        	Xarr=A+B;
        	
        	AA=np.ones([lpd+1, Nz]);
        	AA=self.dz*(0.5+0.5*(self.knots*AA.T));
        	BB=np.ones([lpd+1, Nz]);
        	BB=(BB*self.ZelementArr).T-self.dz;
        	Zarr=AA+BB;
        	
        	self.XArrPlot, self.ZArrPlot=np.meshgrid(Xarr,Zarr);
        	self.VpArrPlot=np.ones( [( (lpd+1) * Nz ), ( (lpd+1) * Nx )] )*Vp;
        	self.VsArrPlot=np.ones( [( (lpd+1) * Nz ), ( (lpd+1) * Nx )] )*Vs;
        	self.RhoArrPlot=np.ones( [( (lpd+1) * Nz ), ( (lpd+1) * Nx )] )*Rho;
        self.regular=True;
        return;
    
    def ReadModel(self, fmodel):
        self.regular=False;
        InArr=np.loadtxt(fmodel);
        self.XArr=InArr[:,0];
        self.ZArr=InArr[:,1];
        self.RhoArr=InArr[:,2];
        self.VpArr=InArr[:,3];
        self.VsArr=InArr[:,4];
        self.Nin=self.RhoArr.size;
        self.xmin=self.XArr.min();
        self.xmax=self.XArr.max();
        self.zmin=self.ZArr.min();
        self.zmax=self.ZArr.max();
        return;
    
    def BlockHomoAnomaly(self, Xmin, Xmax, Zmin, Zmax, Va):
        Xlogic=(self.XArr>=Xmin)*(self.XArr<=Xmax);
        Zlogic=(self.ZArr>=Zmin)*(self.ZArr<=Zmax);
        Logic=Xlogic*Zlogic;
        self.VsArr[Logic]=Va;
        if self.plotflag==True:
            Xlogic=(self.XArrPlot>=Xmin)*(self.XArrPlot<=Xmax);
            Zlogic=(self.ZArrPlot>=Zmin)*(self.ZArrPlot<=Zmax);
            Logic=Xlogic*Zlogic;
            self.VsArrPlot[Logic]=Va;
        return;
    
    def CircleHomoAnomaly(self, Xc, Zc, R, Va):
        dArr = np.sqrt( (self.XArr-Xc)**2 + (self.ZArr-Zc)**2);
        Logic = dArr < R;
        self.VsArr[Logic]=Va;
        if self.plotflag==True:
        	dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2);
        	Logic = dArr < R;
        	self.VsArrPlot[Logic]=Va;
        return;
    
    def CircleGradualAnomaly(self, Xc, Zc, R, Va):
        dArr = np.sqrt( (self.XArr-Xc)**2 + (self.ZArr-Zc)**2);
        dVa = Va - self.Vs;
        delD = R - dArr; 
        self.VsArr = 0.5 * (np.sign(delD) + 1) * delD/R * dVa + self.VsArr;
        if self.plotflag==True:
            dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2);
            dVa = Va - self.Vs;
            delD = R - dArr;
            self.VsArrPlot = 0.5 * (np.sign(delD) + 1) * delD/R * dVa + self.VsArrPlot;
        return;
    
    def WriteModel(self, filename):
        N=self.XArr.size;
        OutArr=np.append(self.XArr, self.ZArr);
        OutArr=np.append(OutArr, self.RhoArr)
        OutArr=np.append(OutArr, self.VpArr);
        OutArr=np.append(OutArr, self.VsArr);
        OutArr=OutArr.reshape( 5, N );
        OutArr=OutArr.T
        np.savetxt(filename, OutArr, fmt='  %.6e')
        return;
    
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
            # knots = np.array([-1.0, -1./7.*np.sqrt(21.), 0.0, 1./7.*np.sqrt(21.), 1.0])
        elif self.lpd == 5:
            knots = np.array([-1.0, -0.7650553239294647, -0.2852315164806451, 0.2852315164806451, 0.7650553239294647, 1.0])
        elif self.lpd == 6:
            knots = np.array([-1.0, -0.8302238962785670, -0.4688487934707142, 0.0, 0.4688487934707142, 0.8302238962785670, 1.0])
        elif self.lpd == 7:
            knots = np.array([-1.0, -0.8717401485096066, -0.5917001814331423,\
                -0.2092992179024789, 0.2092992179024789, 0.5917001814331423, 0.8717401485096066, 1.0])
        self.knots=knots;
        return
    
    def plot(self, ds=1000, unit='km', vmin=3.0, vmax=4.0):
        if self.plotflag==False:
            raise ValueError('No plot array!');
        plt.figure();
        if self.regular==True:
            plt.pcolormesh(self.XArrPlot/ds, self.ZArrPlot/ds, self.VsArrPlot/ds, cmap='seismic_r', vmin=3.0, vmax=4.0);
        else:
            xi = np.linspace(self.xmin, self.xmax, self.Nx*10)
            zi = np.linspace(self.zmin, self.zmax, self.Nz*10)
            self.xi, self.zi = np.meshgrid(xi, zi)
            #-- Interpolating at the points in xi, yi
            self.vi = griddata(self.XArr, self.ZArr, self.VsArr, self.xi, self.zi, 'linear')
            plt.pcolormesh(self.xi/ds, self.zi/ds, ma.getdata(self.vi)/ds, cmap='seismic_r', vmin=3.0, vmax=4.0);
        plt.xlabel(unit);
        plt.ylabel(unit);
        plt.colorbar();
        plt.axis([self.xmin/ds, self.xmax/ds, self.zmin/ds, self.zmax/ds]);
        plt.show()
        return;
    
    def GetMinMaxV(self):
        vmin=self.VsArr.min();
        vmax=self.VsArr.max();
        return vmin, vmax
    
    
class InputChecker(object):
    
    def __init__(self, dt, dx, dz, fc, lpd, vmin, vmax):

        self.dt=dt;
        self.dx=dx;
        self.dz=dz;
        self.fc=fc;
        self.lpd=lpd;
        self.get_GLL();
        self.vmin=vmin;
        self.vmax=vmax;
        return;
    
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
        self.knots=knots;
        return
    
    def CheckMinLambda(self):
        lambdamin=self.vmin/self.fc;
        dxArr=self.dx*np.diff( (0.5+0.5*(self.knots)) );
        dzArr=self.dz*np.diff( (0.5+0.5*(self.knots)) );
        dsmax=max(dxArr.max(), dzArr.max());
        if dsmax * 16. > lambdamin:
            raise ValueError('Grid spacing is too large: ', str(dsmax),' for ',lambdamin, ' km');
        else:
            print 'Grid spacing: ', str(dsmax),'km for ',lambdamin, ' km';
        return
    
    def CheckCFLCondition(self, C=0.35):
        dxArr=self.dx*np.diff( (0.5+0.5*(self.knots)) );
        dzArr=self.dz*np.diff( (0.5+0.5*(self.knots)) );
        dsmin=min(dxArr.min(), dzArr.min());
        dtCFL=C*dsmin/self.vmax;
        
        if self.dt > dtCFL:
            raise ValueError('Time step violates Courant-Frieddrichs-Lewy Condition: ', dt, dtCFL);
        else:
            print 'Time Step: ',self.dt,' s Required Time Step: ', dtCFL, 's ';
        return;
    
    def Check(self, C=0.35):
        print '---------- Checking stability conditions ----------'
        self.CheckMinLambda();
        self.CheckCFLCondition(C=C);
        print '---------- Stability conditions checked  ----------'
        return;
    
    
class WaveSnapshot(object):
    def __init__(self, datadir, nf, ni=5, dn=10, pfx='wavefield', sfx='_01_000.txt',
        gridfname='wavefield_grid_for_dumps_000.txt'):

        self.datadir=datadir;
        self.pfx=pfx;
        self.sfx=sfx;
        self.snapshots=[];
        Ns=int(nf/dn);
        self.Narr=np.arange(Ns)*dn+dn;
        self.Narr=np.append(ni, self.Narr);
        self.gridfname=gridfname;
        return;
    
    def ReadGridFile(self):
        infname=self.datadir+'/'+self.gridfname;
        InArr=np.loadtxt(infname);
        self.xArr=InArr[:,0];
        self.zArr=InArr[:,1];
        self.xmin=self.xArr.min();
        self.xmax=self.xArr.max();
        self.zmin=self.zArr.min();
        self.zmax=self.zArr.max();
        return;
    
    def ReadSnapshots(self):
        wfmax=-999;
        wfmin=999;
        for N in self.Narr:
            infname=(self.datadir+'/'+self.pfx+'%07d'+self.sfx) % (N);
            InArr=np.loadtxt(infname);
            wfmax=max(wfmax, InArr.max() );
            wfmin=min(wfmin, InArr.min() );
            self.wfmax=max(wfmax, abs(wfmin));
            self.snapshots.append(InArr);
        return;
    
    def PlotSnapshots(self, factor=2., outfname=None):
        fig = plt.figure()
        ims = [];
        xi = np.linspace(self.xmin, self.xmax, 100)
        zi = np.linspace(self.zmin, self.zmax, 100)
        self.xi, self.zi = np.meshgrid(xi, zi)
        #-- Interpolating at the points in xi, yi
        i=0;
        for snap in self.snapshots:
            i=i+1
            print 'Plotting ',i,' snapshot!' 
            vi = griddata(self.xArr, self.zArr, snap[::-1], self.xi, self.zi, 'linear')
            im=plt.imshow(vi, cmap='seismic_r', vmin = -self.wfmax/factor, vmax = self.wfmax/factor);
            # im=plt.pcolor(self.xArr, self.zArr, snap)
            ims.append([im])
        
        im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000, blit=True)
        plt.xlabel('x (km)')
        plt.ylabel('z (km)')
        if outfname!=None:
            im_ani.save(outfname, )
        plt.show()
        return im_ani
    
    def SaveDSSnap(self, Nx=100, Nz=100, dpfx='lf'):
        xi = np.linspace(self.xmin, self.xmax, Nx)
        zi = np.linspace(self.zmin, self.zmax, Nz)
        self.xi, self.zi = np.meshgrid(xi, zi)
        #-- Interpolating at the points in xi, yi
        xfname=self.datadir+'/'+dpfx+'_x_'+self.gridfname;
        zfname=self.datadir+'/'+dpfx+'_z_'+self.gridfname;
        np.save(xfname,xi);
        np.save(zfname,zi);
        for i in np.arange(self.Narr.size):
            print 'Downsampling ',i+1,' snapshot!' 
            vi = griddata(self.xArr, self.zArr, self.snapshots[i], self.xi, self.zi, 'linear')
            vfname=(self.datadir+'/'+dpfx+self.pfx+'%07d'+self.sfx) % (self.Narr[i]);
            np.save(vfname,ma.getdata(vi));
        return;
    
    # def SaveDSSnapParallel(self, Nx=100, Nz=100, dpfx='lf'):
    #     xi = np.linspace(self.xmin, self.xmax, Nx)
    #     zi = np.linspace(self.zmin, self.zmax, Nz)
    #     self.xi, self.zi = np.meshgrid(xi, zi)
    #     # -- Interpolating at the points in xi, yi
    #     xfname=self.datadir+'/'+dpfx+'_x_'+self.gridfname;
    #     zfname=self.datadir+'/'+dpfx+'_z_'+self.gridfname;
    #     np.save(xfname,xi);
    #     np.save(zfname,zi);
    #     saveDSsnap = partial(SaveDSsnapshots, xi=self.xi, zi=self.zi, suffix=suffix)
    #     pool =mp.Pool()
    #     pool.map(saveDSsnap, self.snapshots) #make our results with a map call
    #     pool.close() #we are not adding any more processes
    #     pool.join() #tell it to wait until all threads are done before going on
    #     print 'End of Adding Horizontal Slowness  ( Parallel ) !'
    #     
    # 
    #     for i in np.arange(self.Narr.size):
    #         print 'Downsampling ',i,' snapshot!' 
    #         vi = griddata(self.xArr, self.zArr, self.snapshots[i], self.xi, self.zi, 'nn')
    #         vfname=(self.datadir+'/'+dpfx+self.pfx+'%07d'+self.sfx) % (self.Narr[i]);
    #         np.save(vfname,ma.getdata(vi));
    #     return;
    
    def LoadDSSnap(self, dpfx='lf'):
        wfmax=-999;
        wfmin=999;
        #-- Interpolating at the points in xi, yi
        i=0;
        xfname=self.datadir+'/'+dpfx+'_x_'+self.gridfname+'.npy';
        zfname=self.datadir+'/'+dpfx+'_z_'+self.gridfname+'.npy';
        self.xi=np.load(xfname);
        self.zi=np.load(zfname);
        self.xmin=self.xi.min();
        self.xmax=self.xi.max();
        self.zmin=self.zi.min();
        self.zmax=self.zi.max();
        self.dssnaps=[]
        for N in self.Narr:
            i=i+1
            print 'Loading Downsampled ',i,' snapshot!' 
            vfname=(self.datadir+'/'+dpfx+self.pfx+'%07d'+self.sfx+'.npy') % (N);
            snap=np.load(vfname);
            wfmax=max(wfmax, snap.max() );
            wfmin=min(wfmin, snap.min() );
            self.wfmax=max(wfmax, abs(wfmin));
            self.dssnaps.append(snap);
        return;
    
    def PlotDSSnaps(self, ds=1000., factor=4., outfname=None):
        fig = plt.figure(dpi=200)
        ims = [];
        #-- Interpolating at the points in xi, yi
        i=0;
        for snap in self.dssnaps:
            i=i+1
            print 'Plotting ',i,' snapshot!' 
            im=plt.imshow(snap[::-1], cmap='seismic_r', extent=[self.xmin/ds, self.xmax/ds, self.zmin/ds, self.zmax/ds],
                    vmin = -self.wfmax/factor, vmax = self.wfmax/factor);
            # im=plt.pcolor(self.xArr, self.zArr, snap)
            ims.append([im])
        
        im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat=False, blit=True)
        plt.xlabel('x (km)')
        plt.ylabel('z (km)')
        plt.show()
        if outfname!=None:
            try:
                im_ani.save(outfname, dpi=200)
            except:
                return im_ani
        return im_ani



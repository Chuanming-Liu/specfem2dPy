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
        # Generate X grid points given element number and Lagrange polynomial degree
        A=np.ones([(lpd+1)*(lpd+1),Nx]);
        repeatknots=np.tile(self.knots, (lpd+1));
        self.XArr=(A*self.XelementArr).T + (0.5+0.5*repeatknots*A.T)*self.dx - self.dx;
        self.XArr=self.XArr.reshape( (lpd+1)*(lpd+1)*Nx );
        self.XArr=np.tile(self.XArr, Nz)
        # Generate X grid points given element number and Lagrange polynomial degree
        A=np.ones([lpd+1,Nz])
        B=(A*self.ZelementArr).T
        C=(0.5+0.5*A.T*self.knots)*self.dz - self.dz;
        D=B+C
        E=np.tile(D, Nx)
        self.ZArr=np.repeat(E,lpd+1)
        # Generate background velocity and density parameters
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
    
    def plot(self, ds=1000, unit='km', vmin=2.5, vmax=4.5):
        if self.plotflag==False:
            raise ValueError('No plot array!');
        plt.figure(figsize=(16,12));
        if self.regular==True:
            plt.pcolormesh(self.XArrPlot/ds, self.ZArrPlot/ds, self.VsArrPlot/ds, cmap='seismic_r', vmin=vmin, vmax=vmax);
        else:
            xi = np.linspace(self.xmin, self.xmax, self.Nx*10)
            zi = np.linspace(self.zmin, self.zmax, self.Nz*10)
            self.xi, self.zi = np.meshgrid(xi, zi)
            #-- Interpolating at the points in xi, yi
            self.vi = griddata(self.XArr, self.ZArr, self.VsArr, self.xi, self.zi, 'linear')
            plt.pcolormesh(self.xi/ds, self.zi/ds, ma.getdata(self.vi)/ds, cmap='seismic_r', vmin=vmin, vmax=vmax);
        plt.plot( 320, 320 , 'y*', markersize=30)
        plt.xlabel('x('+unit+')', fontsize=30);
        plt.ylabel('z('+unit+')', fontsize=30);
        plt.colorbar();
        plt.axis([self.xmin/ds, self.xmax/ds, self.zmin/ds, self.zmax/ds]);
        # plt.axis('scaled');
        plt.yticks(fontsize=20)
        plt.xticks(fontsize=20)
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
    def __init__(self, datadir, xmax, Nx, zmax, Nz, nf, xmin=0, zmin=0, ni=5, dn=10,
        pfx='wavefield', sfx='_01_000.txt', gridfname='wavefield_grid_for_dumps_000.txt', lpd=4):

        self.datadir=datadir;
        self.pfx=pfx;
        self.sfx=sfx;
        self.snapshots=[];
        Ns=int(nf/dn);
        self.Narr=np.arange(Ns)*dn+dn;
        self.Narr=np.append(ni, self.Narr);
        self.gridfname=gridfname;
        self.dx=(xmax-xmin)/Nx/2;
        self.dz=(zmax-zmin)/Nz/2;
        self.Nx=2*Nx;
        self.Nz=2*Nz;
        self.lpd=lpd;
        XArr = np.arange(2*Nx+1)*self.dx+xmin;
        ZArr = np.arange(2*Nz+1)*self.dz+zmin;
        self.XArr, self.ZArr = np.meshgrid(XArr, ZArr); 
        return;
    
    def ReadGridFile(self):
        print 'Reading Grid File!';
        infname=self.datadir+'/'+self.gridfname;
        InArr=np.loadtxt(infname);
        self.xArrIn=InArr[:,0];
        self.zArrIn=InArr[:,1];
        self.xmin=self.xArrIn.min();
        self.xmax=self.xArrIn.max();
        self.zmin=self.zArrIn.min();
        self.zmax=self.zArrIn.max();
        print 'End Reading Grid File!';
        return;
    
    def GetElementIndex(self):
        print 'Getting element indices !';
        XArr=self.XArr.reshape( (self.Nx+1)*(self.Nz+1) );
        ZArr=self.ZArr.reshape( (self.Nx+1)*(self.Nz+1) );
        self.index=np.array([],dtype=int)
        for i in np.arange( (self.Nx+1)*(self.Nz+1) ):
            x=XArr[i];
            z=ZArr[i];
            Logic = (self.xArrIn==x)*(self.zArrIn==z);
            index=int(np.where(Logic==True)[0][0])
            self.index=np.append(self.index, index)
        print 'End getting element indices !';
        return;
    
    def SaveElementIndex(self, outdir):
        outfname=outdir+'/index.npy';
        np.save(outfname, self.index);
        return;
    
    def LoadElementIndex(self, datadir):
        infname=datadir+'/index.npy';
        self.index=np.load(infname);
        return;
    

    def ReadSnapshots(self):
        wfmax=-999;
        wfmin=999;
        for N in self.Narr:
            infname=(self.datadir+'/'+self.pfx+'%07d'+self.sfx) % (N);
            print 'Reading ',infname,' snapshot!' 
            InArr=np.loadtxt(infname);
            try:
                InArr=InArr[:,1];
            except:
                InArr=InArr;
            wfmax=max(wfmax, InArr.max() );
            wfmin=min(wfmin, InArr.min() );
            self.wfmax=max(wfmax, abs(wfmin));
            snap=np.take(InArr, self.index).reshape(self.Nz+1, self.Nx+1)
            self.snapshots.append( snap[::-1, :] );
        return;
    
    def PlotSnapshots(self, ds=1000., factor=5., outfname=None):
        fig = plt.figure(figsize=(16,12))
        ims = [];
        i=0;
        for snap in self.snapshots:
            i=i+1
            print 'Plotting ',i,' snapshot!' 
            im=plt.imshow(snap, cmap='seismic_r', extent=[self.xmin/ds, self.xmax/ds, self.zmin/ds, self.zmax/ds],
                    vmin = -self.wfmax/factor, vmax = self.wfmax/factor);
            ims.append([im])
        im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000, blit=True)
        plt.xlabel('x (km)')
        plt.ylabel('z (km)')
        if outfname!=None:
            im_ani.save(outfname, )
        plt.show()
        return im_ani
    
    def ReadSingleSnap(self, N):
        wfmax=-999;
        wfmin=999;
        infname=(self.datadir+'/'+self.pfx+'%07d'+self.sfx) % (N);
        print 'Reading ',infname,' snapshot!' 
        InArr=np.loadtxt(infname);
        wfmax=max(wfmax, InArr.max() );
        wfmin=min(wfmin, InArr.min() );
        self.wfmax=max(wfmax, abs(wfmin));
        snap=np.take(InArr, self.index).reshape(self.Nz+1, self.Nx+1)
        self.singleSnap=snap
        return;
    
    def PlotSingleSnap(self, unit='km', ds=1000., factor=1., outfname=None):
        fig = plt.figure(figsize=(16,12))
        im=plt.pcolormesh(self.XArr/ds, self.ZArr/ds, self.singleSnap, shading='gouraud', cmap='seismic_r', vmin = -self.wfmax/factor, vmax = self.wfmax/factor)
        plt.plot( 320, 320 , 'y*', markersize=30)
        plt.xlabel('x('+unit+')', fontsize=30);
        plt.ylabel('z('+unit+')', fontsize=30);
        plt.colorbar();
        plt.axis([self.xmin/ds, self.xmax/ds, self.zmin/ds, self.zmax/ds]);
        # plt.axis('scaled');
        plt.yticks(fontsize=20)
        plt.xticks(fontsize=20)
        plt.show()
        return im
    


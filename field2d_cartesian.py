#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sinter
from matplotlib.mlab import griddata
import numpy.ma as ma
import scipy.ndimage.filters 
from math import pi
# from skimage.filters import roberts, sobel, scharr, prewitt
from scipy.ndimage import convolve
import matplotlib.animation as animation
import pyasdf

X_diff_weight_2 = np.array([[1., 0., -1.]])/2.
Y_diff_weight_2 = X_diff_weight_2.T
X_diff_weight_4 = np.array([[-1., 8., 0., -8., 1.]])/12.
Y_diff_weight_4 = X_diff_weight_4.T
X_diff_weight_6 = np.array([[1./60., 	-3./20.,  3./4.,  0., -3./4., 3./20.,  -1./60.]])
Y_diff_weight_6 = X_diff_weight_6.T

X_diff2_weight_2 = np.array([[1., -2., 1.]])
Y_diff2_weight_2 = X_diff2_weight_2.T
X_diff2_weight_4 = np.array([[-1., 16., -30., 16., -1.]])/12.
Y_diff2_weight_4 = X_diff2_weight_4.T
X_diff2_weight_6 = np.array([[1./90., 	-3./20.,  3./2.,  -49./18., 3./2., -3./20.,  1./90.]])
Y_diff2_weight_6 = X_diff2_weight_6.T


class Field2d(object):
    def __init__(self, Nx, Ny, dx, dy):
        self.Nx=int(Nx)+1
        self.Ny=int(Ny)+1
        self.dx=dx
        self.dy=dy
        self.x=np.arange(Nx+1)*self.dx
        self.y=np.arange(Ny+1)*self.dy
        self.Xarr, self.Yarr = np.meshgrid(self.x, self.y)
        # self.Xarr=self.Xarr.reshape(Nx*Ny)
        # self.Yarr=self.Yarr.reshape(Nx*Ny)
        return
    
    def LoadFile(self, fname):
        try:
            Inarray=np.loadtxt(fname)
        except:
            Inarray=np.load(fname)
        self.XarrIn=Inarray[:,0]/1000.
        self.YarrIn=Inarray[:,1]/1000.
        self.ZarrIn=Inarray[:,2]
        return
    
    def LoadField(self, inField):
        self.XarrIn=inField.Xarr
        self.YarrIn=inField.Yarr
        self.ZarrIn=inField.Zarr
        return
    
    def SaveFile(self, fname, fmt='npy'):
        OutArr=np.append(self.Xarr, self.Yarr)
        OutArr=np.append(OutArr, self.Zarr)
        OutArr=OutArr.reshape(3, self.Nx*self.Ny)
        OutArr=OutArr.T
        if fmt=='npy':
            np.save(fname, OutArr)
        elif fmt=='txt':
            np.savetxt(fname, OutArr)
        else:
            raise TypeError('Wrong output format!')
        return
        
    def interpTest(self, ):
        pass
    
    def Interp(self, kind='cubic', copy=True, bounds_error=False, fill_value=np.nan):
        if self.Xarr.size == self.XarrIn.size:
            if (np.abs(self.Xarr-self.XarrIn)).sum() < 0.01 and (np.abs(self.Yarr-self.YarrIn)).sum() < 0.01:
                print 'No need to interpolate!'
                self.Zarr=self.ZarrIn
                return
        Finterp=sinter.interp2d(self.XarrIn, self.YarrIn, self.ZarrIn,
            kind=kind, copy=copy, bounds_error=bounds_error, fill_value=fill_value)
        self.Zarr = Finterp (self.x, self.y)
        self.Zarr=self.Zarr.reshape(self.Nx*self.Ny)
        return
    
    def natgridInterp(self, interp='linear', copy=True, bounds_error=False, fill_value=np.nan):
        if self.Xarr.size == self.XarrIn.size:
            if (np.abs(self.Xarr.reshape(self.Nx*self.Ny)-self.XarrIn)).sum() < 0.01\
                and (np.abs(self.Yarr.reshape(self.Nx*self.Ny)-self.YarrIn)).sum() < 0.01:
                print 'No need to interpolate!'
                self.Zarr=self.ZarrIn
                return
        # self.Zarr = ma.getdata(griddata(self.XarrIn, self.YarrIn, self.ZarrIn, self.Xarr, self.Yarr, interp=interp))
        self.Zarr = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, self.Xarr, self.Yarr, interp=interp)
        self.Zarr = ma.getdata(self.Zarr)
        # self.Zarr=self.Zarr.reshape(self.Nx*self.Ny)
        return
    
    def PlotInField(self):
        fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        xi, yi = np.meshgrid(self.x, self.y)
        zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(xi, yi, zi, cmap='gist_ncar_r', shading='gouraud')
        levels=np.linspace(zi.min(), zi.max(), 40)
        plt.contour(xi, yi, zi, colors='k', levels=levels)
        plt.axis([0, self.x[-1], 0, self.y[-1]])
        plt.xlabel('km')
        plt.ylabel('km')
        plt.show()
        return
        
    def PlotField(self):
        fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(self.Xarr, self.Yarr, self.Zarr, cmap='gist_ncar_r', shading='gouraud')
        levels=np.linspace(self.Zarr.min(), self.Zarr.max(), 100)
        plt.contour(self.Xarr, self.Yarr, self.Zarr, colors='k', levels=levels)
        plt.axis([0, self.x[-1], 0, self.y[-1]])
        plt.xlabel('km')
        plt.ylabel('km')
        plt.show()
        return
    
    def CuttingEdges(self, nx, ny, fieldtype='TravelT'):
        self.Nx=self.Nx-2*nx
        self.Ny=self.Ny-2*ny
        self.x=np.arange(self.Nx)*self.dx
        self.y=np.arange(self.Ny)*self.dy
        self.Xarr, self.Yarr = np.meshgrid(self.x, self.y)
        self.Zarr=self.Zarr[nx:-nx, ny:-ny]
        return
        
    def fftDiff(self, m, n):
        try:
            h = np.fft.fft2(ma.getdata(self.Zarr))
        except:
            h = np.fft.fft2(self.Zarr)
        hshift = np.fft.fftshift(h)
        Nx=self.Nx
        Ny=self.Ny
        if Nx % 2 ==0:
            u=np.arange(Nx) - Nx/2.
        else:
            u=np.arange(Nx) - (Nx-1)/2.
        if Ny % 2 ==0:
            v=np.arange(Ny) - Ny/2.
        else:
            v=np.arange(Ny) - (Ny-1)/2.
        U,V=np.meshgrid(u,v)
        hdiff =  ((1j*2*np.pi*U/Nx)**m)*((1j*2*np.pi*V/Ny)**n) * hshift
        out_diff = np.real( np.fft.ifft2( np.fft.ifftshift(hdiff) ) )/(self.dx**m)/(self.dy**n)
        out_diff = out_diff[:self.Ny, :self.Nx]
        return out_diff
    
    def fftDiff2(self, m, n):
        Nx=1<<(self.Nx-1).bit_length()
        Ny=1<<(self.Ny-1).bit_length()
        # h = np.fft.fft2(self.Zarr, s=[Nx, Ny] )
        h = np.fft.fft2(ma.getdata(self.Zarr), s=[Nx, Ny] )
        hshift = np.fft.fftshift(h)
        u = np.arange(Nx) - Nx/2.
        v = np.arange(Ny) - Ny/2.
        U,V = np.meshgrid(u,v)
        hdiff = ( (1j*2*np.pi*U/Nx)**m )*( (1j*2*np.pi*V/Ny)**n ) * hshift
        out_diff = np.real( np.fft.ifft2( np.fft.ifftshift(hdiff) ) )/(self.dx**m)/(self.dy**n)
        out_diff = out_diff[:self.Ny, :self.Nx]
        return out_diff
    
    
    def Gradient(self, edge_order=1, method='default', order=2):
        if method=='default':
            self.grad=np.gradient( ma.getdata(self.Zarr), self.dx, self.dy, edge_order=edge_order)
        elif method=='freq':
            diff_x=self.fftDiff(m=1, n=0)
            diff_y=self.fftDiff(m=0, n=1)
            self.grad=[]
            self.grad.append(diff_y)
            self.grad.append(diff_x)
        elif method == 'convolve':
            if order==2:
                diff_x=convolve( ma.getdata(self.Zarr), X_diff_weight_2)/self.dx
                diff_y=convolve(ma.getdata(self.Zarr), Y_diff_weight_2)/self.dy
            elif order==4:
                diff_x=convolve(ma.getdata(self.Zarr), X_diff_weight_4)/self.dx
                diff_y=convolve(ma.getdata(self.Zarr), Y_diff_weight_4)/self.dy
            elif order==6:
                diff_x=convolve(ma.getdata(self.Zarr), X_diff_weight_6)/self.dx
                diff_y=convolve(ma.getdata(self.Zarr), Y_diff_weight_6)/self.dy
            self.grad=[]
            self.grad.append(diff_y)
            self.grad.append(diff_x)
        return
    
    def GetApparentV(self):
        self.AppV = np.sqrt ( self.grad[0] ** 2 + self.grad[1] ** 2)
        self.AppV[ np.where(self.AppV==0) ] = -1.
        self.AppV=1./self.AppV
        return
    
    def PlotAppV(self, vmin=2.9, vmax=3.1):
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        plt.subplots()
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        # plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.AppV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax)
        plt.pcolormesh(self.Xarr, self.Yarr, self.AppV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax)
        # plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.AppV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax)
        plt.axis('equal')
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]])
        plt.xlabel('km')
        plt.ylabel('km')
        # plt.show()
        return
    
    def GetDoT(self, enx, eny):
        Darr=np.sqrt( (self.Xarr-enx)**2 + (self.Yarr-eny)**2)
        self.DoT = Darr/ma.getdata(self.Zarr)
        return
    
    def PlotDoT(self, vmin=2.9, vmax=3.1):
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        plt.subplots()
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        # plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.AppV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax)
        plt.pcolormesh(self.Xarr, self.Yarr, self.DoT, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax)
        plt.axis('equal')
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]])
        plt.xlabel('km')
        plt.ylabel('km')
        # plt.show()
        return
    
    def LaplacianEqualXY(self):
        if self.dx!=self.dy:
            raise ValueError('grid spacing not equal!')
        self.lplc=scipy.ndimage.filters.laplace(ma.getdata(self.Zarr) ) / (self.dx*self.dy)
        self.lplc=self.lplc[1:-1, 1:-1]
        return
    
    def Laplacian(self, method='default', order=2):
        Zarr=ma.getdata(self.Zarr)
        if method == 'default':
            Zarr_yp=Zarr[2:, 1:-1]
            Zarr_yn=Zarr[:-2, 1:-1]
            Zarr_xp=Zarr[1:-1, 2:]
            Zarr_xn=Zarr[1:-1, :-2]
            Zarr=Zarr[1:-1, 1:-1]
            self.lplc=(Zarr_yp+Zarr_yn-2*Zarr) / (self.dy**2) + (Zarr_xp+Zarr_xn-2*Zarr) / (self.dx**2)
        elif method == 'convolve':
            if order==2:
                diff2_x=convolve( ma.getdata(self.Zarr), X_diff2_weight_2)/self.dx/self.dx
                diff2_y=convolve(ma.getdata(self.Zarr), Y_diff2_weight_2)/self.dy/self.dy
            elif order==4:
                diff2_x=convolve( ma.getdata(self.Zarr), X_diff2_weight_4)/self.dx/self.dx
                diff2_y=convolve(ma.getdata(self.Zarr), Y_diff2_weight_4)/self.dy/self.dy
            elif order==6:
                diff2_x=convolve( ma.getdata(self.Zarr), X_diff2_weight_6)/self.dx/self.dx
                diff2_y=convolve(ma.getdata(self.Zarr), Y_diff2_weight_6)/self.dy/self.dy
            self.lplc=diff2_x+diff2_y
            self.lplc=self.lplc[1:-1, 1:-1]
        return
    
    
    
    def PlotLaplacian(self):
        fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.lplc, cmap='seismic', shading='gouraud' )
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]])
        plt.xlabel('km')
        plt.ylabel('km')
        # plt.show()
        return
    
    def GetLplcCorrection(self, per):
        omega=2.*np.pi/per
        Zarr=ma.getdata(self.Zarr)
        self.lplcCo=self.lplc/Zarr[1:-1, 1:-1]/(omega**2)
        return
    
    def PlotLplcCo(self):
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        plt.subplots()
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.lplcCo, cmap='seismic', shading='gouraud' , vmin=-0.01, vmax=0.01)
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]])
        plt.xlabel('km')
        plt.ylabel('km')
        # plt.show()
        return
    
    def GetCorV(self, inAmpField):
        lplcCo=inAmpField.lplcCo
        try:
            self.CorV=1./ np.sqrt(1./(self.AppV**2) - lplcCo)
        except:
            self.CorV=1./ np.sqrt(1./(self.AppV[1:-1, 1:-1]**2) - lplcCo)
        return
    
    def PlotCorV(self, vmin=2.9, vmax=3.1):
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        plt.subplots()
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.CorV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax)
        plt.axis('equal')
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]])
        plt.xlabel('km')
        plt.ylabel('km')
        # plt.show()
        return
    
    
class WaveSnapshot(object):
    """
    An object to handle output wavefield from SPECFEM2D.
    --------------------------------------------------------------------------------------------------------------
    Parameters:
    datadir            - data directory
    pfx                   - input file prefix
    sfx                    - input file suffix
    gridfname       - grid point file name
    snapshots        - snapshot list
    Narr                - snapshot number array
    dx, dz              - (half) element size
    Nx, Nz            - (doubled) element number in x, z 
                                note that some data points locate at the mid point between two element point 
    XArr, ZArr    - arrays for element location
    lpd                   - Lagrange polynomial degree
    --------------------------------------------------------------------------------------------------------------
    """
    def __init__(self, datadir, xmax, Nx, zmax, Nz, nt, xmin=0, zmin=0, ni=5, dn=10,
        pfx='wavefield', sfx='_01_000.txt', gridfname='wavefield_grid_for_dumps_000.txt', lpd=4):
        """ 
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        xmin, zmax    - bound of study region
        nt                     - number of total time step
        ni                     - initial time step number
        dn                    - time step interval
        -----------------------------------------------------------------------------------------------------
        """
        self.datadir=datadir
        self.pfx=pfx
        self.sfx=sfx
        self.snapshots=[]
        Ns=int(nt/dn)
        self.Narr=np.arange(Ns)*dn+dn
        self.Narr=np.append(ni, self.Narr)
        self.gridfname=gridfname
        self.dx=(xmax-xmin)/Nx/2 #  half of element size
        self.dz=(zmax-zmin)/Nz/2
        self.Nx=2*Nx # double the element number
        self.Nz=2*Nz
        self.lpd=lpd
        XArr = np.arange(2*Nx+1)*self.dx+xmin
        ZArr = np.arange(2*Nz+1)*self.dz+zmin
        self.XArr, self.ZArr = np.meshgrid(XArr, ZArr) 
        return
    
    def ReadGridFile(self):
        """
        Read grid point file
        """
        print 'Reading Grid File !'
        infname=self.datadir+'/'+self.gridfname
        InArr=np.loadtxt(infname)
        self.xArrIn=InArr[:,0]
        self.zArrIn=InArr[:,1]
        self.xmin=self.xArrIn.min()
        self.xmax=self.xArrIn.max()
        self.zmin=self.zArrIn.min()
        self.zmax=self.zArrIn.max()
        print 'End Reading Grid File !'
        return
    
    def GetElementIndex(self):
        """
        Get the element indices
        """
        print 'Getting element indices !'
        Ntotal=(self.Nx+1)*(self.Nz+1)
        XArr=self.XArr.reshape( Ntotal )
        ZArr=self.ZArr.reshape( Ntotal )
        self.index=np.array([],dtype=int)
        for i in np.arange( Ntotal ):
            if i%1000==0:
                print 'Step:', i, 'of', Ntotal
            x=XArr[i]
            z=ZArr[i]
            Logic = (self.xArrIn==x)*(self.zArrIn==z)
            index=int( np.where( Logic==True)[0][0] )
            self.index=np.append(self.index, index)
        print 'End getting element indices !'
        return
    
    def SaveElementIndex(self, outdir):
        """
        Save the element indices
        """
        outfname=outdir+'/index.npy'
        np.save(outfname, self.index)
        return
    
    def LoadElementIndex(self, datadir):
        """
        Load the element indices
        """
        infname=datadir+'/index.npy'
        self.index=np.load(infname)
        return
    
    def ReadSnapshots(self):
        """
        Read snapshots
        """
        wfmax=-999
        wfmin=999
        for N in self.Narr:
            infname=(self.datadir+'/'+self.pfx+'%07d'+self.sfx) % (N)
            print 'Reading ',infname,' snapshot!' 
            InArr=np.loadtxt(infname)
            try:
                InArr=InArr[:,1]
            except:
                InArr=InArr
            wfmax=max(wfmax, InArr.max() )
            wfmin=min(wfmin, InArr.min() )
            self.wfmax=max(wfmax, abs(wfmin))
            snap=np.take(InArr, self.index).reshape(self.Nz+1, self.Nx+1)
            self.snapshots.append( snap[::-1, :] )
        return
    
    def writeASDF(self, outfname):
        """
        Write wavefield snapshots to ASDF Dataset
        """
        dbase=pyasdf.ASDFDataSet(outfname)
        for i in np.arange(self.Narr.size):
            snap = self.snapshots[i]
            path='wf'+str(int(self.Narr[i]))
            # print path
            wfmax=max(snap.max(), abs(snap.min()))
            header={'xmin': self.xmin, 'xmax': self.xmax, 'zmin': self.zmin, 'zmax': self.zmax, 'wfmaxglobe': self.wfmax, 'wfmax': wfmax}
            dbase.add_auxiliary_data(data=snap, data_type='Snapshot', path=path, parameters=header)
        return
            
    def readASDF(self, infname):
        """
        Read wavefield snapshots from ASDF Dataset
        """
        dbase=pyasdf.ASDFDataSet(infname)
        for i in np.arange(self.Narr.size):
            path='wf'+str(int(self.Narr[i]))
            snap=dbase.auxiliary_data.Snapshot[path].data.value
            self.snapshots.append(snap)
        self.wfmax=dbase.auxiliary_data.Snapshot[path].parameters['wfmaxglobe']
        self.xmin=dbase.auxiliary_data.Snapshot[path].parameters['xmin']
        self.xmax=dbase.auxiliary_data.Snapshot[path].parameters['xmax']
        self.zmin=dbase.auxiliary_data.Snapshot[path].parameters['zmin']
        self.zmax=dbase.auxiliary_data.Snapshot[path].parameters['zmax']
        return
    
    
    def PlotSnapshots(self, ds=1000., factor=5., outfname=None):
        """
        Plot snapshots as animation
        """
        fig = plt.figure(figsize=(16,12))
        ims = []
        i=0
        for snap in self.snapshots:
            i=i+1
            print 'Plotting ',i,' snapshot!' 
            im=plt.imshow(snap, cmap='seismic_r', extent=[self.xmin/ds, self.xmax/ds, self.zmin/ds, self.zmax/ds],
                    vmin = -self.wfmax/factor, vmax = self.wfmax/factor)
            ims.append([im])
        im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000, blit=True)
        plt.xlabel('x (km)')
        plt.ylabel('z (km)')
        if outfname!=None:
            im_ani.save(outfname, )
        plt.show()
        return im_ani
    
    def ReadSingleSnap(self, N):
        """
        Read single snapshot
        """
        wfmax=-999
        wfmin=999
        infname=(self.datadir+'/'+self.pfx+'%07d'+self.sfx) % (N)
        print 'Reading ',infname,' snapshot!' 
        InArr=np.loadtxt(infname)
        wfmax=max(wfmax, InArr.max() )
        wfmin=min(wfmin, InArr.min() )
        self.wfmax=max(wfmax, abs(wfmin))
        snap=np.take(InArr, self.index).reshape(self.Nz+1, self.Nx+1)
        self.singleSnap=snap
        return
    
    def PlotSingleSnap(self, unit='km', ds=1000., factor=1., outfname=None):
        """
        Plot single snapshot
        """
        fig = plt.figure(figsize=(16,12))
        im=plt.pcolormesh(self.XArr/ds, self.ZArr/ds, self.singleSnap, shading='gouraud', cmap='seismic_r', vmin = -self.wfmax/factor, vmax = self.wfmax/factor)
        plt.plot( 320, 320 , 'y*', markersize=30)
        plt.xlabel('x('+unit+')', fontsize=30)
        plt.ylabel('z('+unit+')', fontsize=30)
        plt.colorbar()
        plt.axis([self.xmin/ds, self.xmax/ds, self.zmin/ds, self.zmax/ds])
        # plt.axis('scaled')
        plt.yticks(fontsize=20)
        plt.xticks(fontsize=20)
        plt.show()
        return im
    


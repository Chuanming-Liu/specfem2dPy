#!/usr/bin/env python
import numpy as np;
import matplotlib.pyplot as plt;
import scipy.interpolate as sinter
from matplotlib.mlab import griddata
import numpy.ma as ma
import scipy.ndimage.filters 


class Field2d(object):
    def __init__(self, Nx, Ny, dx, dy):
        self.Nx=int(Nx)+1;
        self.Ny=int(Ny)+1;
        self.dx=dx;
        self.dy=dy;
        self.x=np.arange(Nx+1)*self.dx;
        self.y=np.arange(Ny+1)*self.dy;
        self.Xarr, self.Yarr = np.meshgrid(self.x, self.y)
        # self.Xarr=self.Xarr.reshape(Nx*Ny);
        # self.Yarr=self.Yarr.reshape(Nx*Ny);
        return;
    
    def LoadFile(self, fname):
        try:
            Inarray=np.loadtxt(fname);
        except:
            Inarray=np.load(fname);
        self.XarrIn=Inarray[:,0]/1000.;
        self.YarrIn=Inarray[:,1]/1000.;
        self.ZarrIn=Inarray[:,2];
        return;
    
    def LoadField(self, inField):
        self.XarrIn=inField.Xarr;
        self.YarrIn=inField.Yarr;
        self.ZarrIn=inField.Zarr;
        return;
    
    def SaveFile(self, fname, fmt='npy'):
        OutArr=np.append(self.Xarr, self.Yarr);
        OutArr=np.append(OutArr, self.Zarr);
        OutArr=OutArr.reshape(3, self.Nx*self.Ny);
        OutArr=OutArr.T;
        if fmt=='npy':
            np.save(fname, OutArr);
        elif fmt=='txt':
            np.savetxt(fname, OutArr);
        else:
            raise TypeError('Wrong output format!')
        return;
        
    
    def Interp(self, kind='cubic', copy=True, bounds_error=False, fill_value=np.nan):
        if self.Xarr.size == self.XarrIn.size:
            if (np.abs(self.Xarr-self.XarrIn)).sum() < 0.01 and (np.abs(self.Yarr-self.YarrIn)).sum() < 0.01:
                print 'No need to interpolate!'
                self.Zarr=self.ZarrIn;
                return;
        Finterp=sinter.interp2d(self.XarrIn, self.YarrIn, self.ZarrIn,
            kind=kind, copy=copy, bounds_error=bounds_error, fill_value=fill_value);
        self.Zarr = Finterp (self.x, self.y);
        self.Zarr=self.Zarr.reshape(self.Nx*self.Ny);
        return;
    
    def natgridInterp(self, interp='linear', copy=True, bounds_error=False, fill_value=np.nan):
        if self.Xarr.size == self.XarrIn.size:
            if (np.abs(self.Xarr.reshape(self.Nx*self.Ny)-self.XarrIn)).sum() < 0.01\
                and (np.abs(self.Yarr.reshape(self.Nx*self.Ny)-self.YarrIn)).sum() < 0.01:
                print 'No need to interpolate!'
                self.Zarr=self.ZarrIn;
                return;
        # self.Zarr = ma.getdata(griddata(self.XarrIn, self.YarrIn, self.ZarrIn, self.Xarr, self.Yarr, interp=interp));
        self.Zarr = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, self.Xarr, self.Yarr, interp=interp);
        # self.Zarr=self.Zarr.reshape(self.Nx*self.Ny);
        return;
    
    def PlotInField(self):
        fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k');
        xi, yi = np.meshgrid(self.x, self.y)
        zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(xi, yi, zi, cmap='gist_ncar_r', shading='gouraud');
        levels=np.linspace(zi.min(), zi.max(), 40)
        plt.contour(xi, yi, zi, colors='k', levels=levels)
        plt.axis([0, self.x[-1], 0, self.y[-1]]);
        plt.xlabel('km');
        plt.ylabel('km');
        plt.show();
        return;
        
    def PlotField(self):
        fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k');
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(self.Xarr, self.Yarr, self.Zarr, cmap='gist_ncar_r', shading='gouraud');
        levels=np.linspace(self.Zarr.min(), self.Zarr.max(), 100)
        plt.contour(self.Xarr, self.Yarr, self.Zarr, colors='k', levels=levels)
        plt.axis([0, self.x[-1], 0, self.y[-1]]);
        plt.xlabel('km');
        plt.ylabel('km');
        plt.show();
        return;
        
    def Gradient(self, edge_order=1):
        # Zarr=ma.getdata(self.Zarr);
        # Zarr_yp=Zarr[2:, 1:-1];
        # Zarr_yn=Zarr[:-2, 1:-1];
        # Zarr_xp=Zarr[1:-1, 2:];
        # Zarr_xn=Zarr[1:-1, :-2];
        # Zarr=Zarr[1:-1, 1:-1];
        # self.grad=[];
        # self.grad.append ( (Zarr_xp-Zarr_xn) / (2.*self.dx) )
        # self.grad.append ( (Zarr_yp-Zarr_yn) / (2.*self.dy) )
        
        self.grad=np.gradient( ma.getdata(self.Zarr), self.dx, self.dy, edge_order=edge_order)
        # self.grad[0]=self.grad[0][1:-1, 1:-1];
        # self.grad[1]=self.grad[1][1:-1, 1:-1];
        return;
    
    def GetApparentV(self):
        self.AppV = np.sqrt ( self.grad[0] ** 2 + self.grad[1] ** 2);
        self.AppV[ np.where(self.AppV==0) ] = -1.;
        self.AppV=1./self.AppV;
        return;
    
    def PlotAppV(self, vmin=3.4, vmax=3.6):
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k');
        plt.subplots()
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        # plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.AppV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax);
        plt.pcolormesh(self.Xarr, self.Yarr, self.AppV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax);
        # plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.AppV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax);
        plt.axis('equal')
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]]);
        plt.xlabel('km');
        plt.ylabel('km');
        # plt.show();
        return;
    
    def GetDoT(self):
        Darr=np.sqrt( (self.Xarr-240)**2 + (self.Yarr-240)**2);
        self.DoT = Darr/ma.getdata(self.Zarr);
        return;
    
    def PlotDoT(self, vmin=3.4, vmax=3.6):
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k');
        plt.subplots()
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        # plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.AppV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax);
        plt.pcolormesh(self.Xarr, self.Yarr, self.DoT, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax);
        plt.axis('equal')
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]]);
        plt.xlabel('km');
        plt.ylabel('km');
        # plt.show();
        return;
    
    def LaplacianEqualXY(self):
        if self.dx!=self.dy:
            raise ValueError('grid spacing not equal!');
        self.lplc=scipy.ndimage.filters.laplace(ma.getdata(self.Zarr) ) / (self.dx*self.dy);
        self.lplc=self.lplc[1:-1, 1:-1];
        return;
    
    def Laplacian(self):
        Zarr=ma.getdata(self.Zarr);
        Zarr_yp=Zarr[2:, 1:-1];
        Zarr_yn=Zarr[:-2, 1:-1];
        Zarr_xp=Zarr[1:-1, 2:];
        Zarr_xn=Zarr[1:-1, :-2];
        Zarr=Zarr[1:-1, 1:-1]
        self.lplc=(Zarr_yp+Zarr_yn-2*Zarr) / (self.dy**2) + (Zarr_xp+Zarr_xn-2*Zarr) / (self.dx**2);
        return;
    
    def PlotLaplacian(self):
        fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k');
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.lplc, cmap='seismic', shading='gouraud' );
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]]);
        plt.xlabel('km');
        plt.ylabel('km');
        plt.show();
        return;
    
    def GetLplcCorrection(self, per):
        omega=2.*np.pi/per;
        Zarr=ma.getdata(self.Zarr)
        self.lplcCo=self.lplc/Zarr[1:-1, 1:-1]/(omega**2);
        return
    
    def PlotLplcCo(self):
        fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k');
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.lplcCo, cmap='seismic', shading='gouraud' , vmin=-0.01, vmax=0.01);
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]]);
        plt.xlabel('km');
        plt.ylabel('km');
        # plt.show();
        return;
    
    def GetCorV(self, inAmpField):
        lplcCo=inAmpField.lplcCo;
        self.CorV=1./ np.sqrt(1./(self.AppV**2) - lplcCo);
        return;
    
    def PlotCorV(self, vmin=1.5, vmax=4.5):
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k');
        plt.subplots()
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.CorV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax);
        plt.axis('equal')
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]]);
        plt.xlabel('km');
        plt.ylabel('km');
        # plt.show();
        return;
    
    
    
    
    
    # def Laplacian(self, ):
        
    
    
    #         cb=m.colorbar(location='bottom',size='2%')
    # cb.ax.tick_params(labelsize=15)
    # cb.set_label('Phase travel time (sec)', fontsize=20)
    
        
        
            
        
        
    
        
    
    


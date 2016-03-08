#!/usr/bin/env python
"""
2D wave equation solver based on fatiando.
"""
import obspy
from fatiando.seismic import wavefd
import numpy as np
from fatiando.vis import mpl
import matplotlib.pyplot as plt
from matplotlib import animation
import os

#-------------------------
def make_colormap(colors):
#-------------------------
    """
    Define a new color map based on values specified in the dictionary
    colors, where colors[z] is the color that value z should be mapped to,
    with linear interpolation between the given values of z.
    The z values (dictionary keys) are real numbers and the values
    colors[z] can be either an RGB list, e.g. [1,0,0] for red, or an
    html hex string, e.g. "#ff0000" for red.
    """
    from matplotlib.colors import LinearSegmentedColormap, ColorConverter
    from numpy import sort
    z = sort(colors.keys())
    n = len(z)
    z1 = min(z)
    zn = max(z)
    x0 = (z - z1) / (zn - z1)
    CC = ColorConverter()
    R = []
    G = []
    B = []
    for i in range(n):
        #i'th color at level z[i]:
        Ci = colors[z[i]]
        if type(Ci) == str:
            # a hex string of form '#ff0000' for example (for red)
            RGB = CC.to_rgb(Ci)
        else:
            # assume it's an RGB triple already:
            RGB = Ci
        R.append(RGB[0])
        G.append(RGB[1])
        B.append(RGB[2])
    cmap_dict = {}
    cmap_dict['red'] = [(x0[i],R[i],R[i]) for i in range(len(R))]
    cmap_dict['green'] = [(x0[i],G[i],G[i]) for i in range(len(G))]
    cmap_dict['blue'] = [(x0[i],B[i],B[i]) for i in range(len(B))]
    mymap = LinearSegmentedColormap('mymap',cmap_dict)
    return mymap

class Receivers(object):
    
    def __init__(self, shape, ds=100.):
        self.shape=shape;
        self.ds=ds;
        self.StaLst=[];
        return;
    
    def SingleReceiver(self, Nx, Nz):
        self.StaLst.append([Nx*self.ds, Nz*self.ds]);
    
    def HomoStaLst(self, dx=1, dz=1):
        Nx = int(self.shape[0]/dx);
        Nz = int(self.shape[1]/dz);
        # print Nx, Nz
        for i in np.arange(Nx):
            for j in np.arange(Nz):
                X=i*dx*self.ds;
                Z=j*dz*self.ds;
                newSta=[X, Z];
                self.StaLst.append(newSta);
        return;
    
class Sources(object):
    
    def __init__(self, shape, ds=100.):
        self.shape=shape;
        self.ds=ds;
        self.area = [0, (shape[0]-1) * ds, 0, (shape[1]-1) * ds ];
        self.EventLst=[];
        self.STFLst=[];
        return;
    
    def SingleSource(self, Nx, Nz, fc, amp=1. ):
        self.STFLst.append(wavefd.GaussSourceLF(Nx * self.ds, Nz * self.ds, self.area, self.shape,  amp, fc))
        self.EventLst.append([Nx * self.ds, Nz * self.ds])
        return
    # def HomoStaLst(self, dx=1, dz=1):
    #     Nx = int(shape[0]/dx);
    #     Nz = int(shape[1]/dz);
    #     for i in np.arange(Nx):
    #         for j in np.arange(Nz):
    #             X=i*dx*self.ds;
    #             Z=j*dz*self.ds;
    #             newSta=[X, Z];
    #             self.StaLst.append(newSta);
    #     return;

class VelocityModel(object):
    
    def __init__(self, shape, ds, Vback=3000.):

        self.shape=shape;
        self.ds=ds;
        self.velocity = np.zeros(shape) + Vback;
        self.Vback=Vback
        return;
    
    def BlockHomoAnomaly(self, Xmin, Xmax, Zmin, Zmax, Va):
        if Xmax > self.shape[0] or Zmax > self.shape[1]:
            raise ValueError('Maximum x/z node number exceed study region limit!');
        self.velocity[Xmin:Xmax, Zmin:Zmax] = Va;
        return;
    
    def CircleHomoAnomaly(self, Xc, Zc, R, Va):
        Xmin = Xc - R;
        Xmax = Xc + R;
        Zmin = Zc - R;
        Zmax = Zc + R;
        if Xmax > self.shape[0] or Zmax > self.shape[1] or Xmin < 0 or Zmin < 0:
            raise ValueError('Maximum x/z node number exceed study region limit!');
        for i in np.arange(2*R)+Xmin:
            for j in np.arange(2*R)+Zmin:
                if (i-Xc)**2+(j-Zc)**2 > R*R:
                    continue;
                self.velocity[i, j] = Va;
        return;
    
    def CircleGradualAnomaly(self, Xc, Zc, R, Va):
        Xmin = Xc - R;
        Xmax = Xc + R;
        Zmin = Zc - R;
        Zmax = Zc + R;
        if Xmax > self.shape[0] or Zmax > self.shape[1] or Xmin < 0 or Zmin < 0:
            raise ValueError('Maximum x/z node number exceed study region limit!');
        Vref=self.velocity[Xmin, Zmin];
        for i in np.arange(2*R)+Xmin:
            for j in np.arange(2*R)+Zmin:
                myR = np.sqrt( (i-Xc)**2+(j-Zc)**2 );
                if myR > R:
                    continue;
                self.velocity[i, j] = Va+(Vref-Va)*(myR/R);
        return;
    
    def WriteModel(self, filename):
        np.save(filename, self.velocity);
        return;
    
    def ReadModel(self, filename):
        self.velocity=np.load(filename);
        if self.velocity.shape != self.shape:
            self.shape = self.velocity.shape;
            print 'WARNING: Loaded velocity model has different shape!';
        return;
    
    
class InputDataBase(object):
    
    def __init__(self, shape, ds, dt, duration, padding, sources, fmin, receivers, velocityM, snapshots=None, PlotRec=[]):
        
        self.shape = shape;
        self.ds = ds;
        self.area = [0, (shape[0]-1) * ds, 0, (shape[1]-1) * ds];
        maxdt = wavefd.scalar_maxdt(self.area, self.shape, np.max(velocityM.velocity));
        if dt > maxdt:
            print 'WARNING: Time step: ' +str(dt)+'s violates Courant Condition(' + str(maxdt)+'s )!'
            dt = maxdt;
        self.dt = dt;
        self.duration = duration;
        self.maxit = int(duration / dt);
        self.snapshots = snapshots;
        self.padding = padding;
        self.receivers = receivers;
        self.velocityM = velocityM; 
        self.sources = sources;
        self.PlotRec=PlotRec;
        lambdamin = ( velocityM.velocity.min() ) / fmin ;
        if ds > lambdamin / 16:
            raise ValueError('Required spacing is : '+str(lambdamin / 16)+ ' m ! ');
        return;



class WaveSolverMembrane(object):
    
    def __init__(self, INDataBase):
        self.INDataBase = INDataBase;
        return;
    
    def Run(self, vmin=-20000, vmax=20000):
        INDataBase = self.INDataBase;
        if INDataBase.snapshots is None:
            self.simulation = wavefd.membrane(INDataBase.velocityM.velocity, INDataBase.area, INDataBase.dt, INDataBase.maxit,
                INDataBase.sources.STFLst, INDataBase.receivers.StaLst, INDataBase.snapshots, INDataBase.padding);
            self.t, self.wavefield, self.seismograms = self.simulation.next();
            print 'ATTENTION: End of Computing!'
        else:
            global simulation, wavefield;
            simulation = wavefd.membrane(INDataBase.velocityM.velocity, INDataBase.area, INDataBase.dt, INDataBase.maxit,
                INDataBase.sources.STFLst, INDataBase.receivers.StaLst, INDataBase.snapshots, INDataBase.padding);
            # This part makes an animation using matplotlibs animation API
            background = (INDataBase.velocityM.velocity - INDataBase.velocityM.Vback)/ INDataBase.velocityM.Vback
            fig = mpl.figure(figsize=(16, 12))
            if len(INDataBase.PlotRec) != 0:
                Nplot=len(INDataBase.PlotRec)+1
                ax=mpl.subplot(Nplot,1,1)
                # ax = plt.subplot2grid((Nplot,1),(0, 2), rowspan=2)
                # ax=mpl.subplots()
                ax.xaxis.tick_top()
                ax.xaxis.set_label_position('top') 
                for i in INDataBase.PlotRec:
                    mpl.points(INDataBase.receivers.StaLst[i], '^b', size=8);
            # mpl.subplots_adjust(right=0.98, left=0.11, hspace=0.5, top=0.93)
            if not (background.max()<0.0001 and background.min() > -0.0001):
                Vbackground = plt.contour(background, extent=INDataBase.area,
                                       cmap=mpl.cm.seismic_r, vmin=-1, vmax=1)
            wavefield= mpl.imshow(np.zeros_like(INDataBase.velocityM.velocity), extent=INDataBase.area,
                                   cmap=mpl.cm.seismic_r, vmin=vmin, vmax=vmax);
            # mpl.points(INDataBase.receivers.StaLst, '^b', size=8);
            mpl.points(INDataBase.sources.EventLst, '*g', size=8)
            mpl.ylim(INDataBase.area[2:][::-1])
            mpl.xlabel('x (km)')
            mpl.ylabel('z (km)')
            mpl.m2km()
            if len(INDataBase.PlotRec) != 0:
                # ax=mpl.subplots()
                global Plot_seismograms
                global seis_index
                seis_index=INDataBase.PlotRec;
                Plot_seismograms=[];
                global times
                times = np.linspace(0, INDataBase.dt * INDataBase.maxit, INDataBase.maxit)
                for i in range(len(INDataBase.PlotRec)):
                    mpl.subplot(Nplot,1,i+2)
                    # plt.subplot2grid((Nplot,1),(0, 2), rowspan=2)
                    seismogram1, = mpl.plot([], [], '-k', lw=3);
                    Plot_seismograms.append(seismogram1)
                    mpl.xlim(0, INDataBase.duration)
                    mpl.ylim(vmin/2, vmax/2)
                    mpl.title('X: '+str(INDataBase.receivers.StaLst[INDataBase.PlotRec[i]][0]/1000.)+'km Z:'+
                              str(INDataBase.receivers.StaLst[INDataBase.PlotRec[i]][1]/1000.)+'km')
                ani = animation.FuncAnimation( fig, animate1,  frames=INDataBase.maxit / INDataBase.snapshots, interval=1)
            else:
                ani = animation.FuncAnimation( fig, animate2,  interval=1)
            mpl.show()
            # self.seismograms = seismograms;
        return
    
    def Run_sh(self, density=2500., vmin=-10 ** (-5), vmax=10 ** (-5), outmp4=''):
        INDataBase = self.INDataBase;
        denArr=np.ones(INDataBase.velocityM.velocity.shape)*density;
        muArr=wavefd.lame_mu(INDataBase.velocityM.velocity, denArr)
        if INDataBase.snapshots is None:
            self.simulation = wavefd.membrane_sh(muArr, denArr, INDataBase.area, INDataBase.dt, INDataBase.maxit,
                INDataBase.sources.STFLst, INDataBase.receivers.StaLst, INDataBase.snapshots, INDataBase.padding);
            self.t, self.wavefield, self.seismograms = self.simulation.next();
            print 'ATTENTION: End of Computing!'
        else:
            global simulation, wavefield;
            simulation = wavefd.membrane_sh(muArr, denArr, INDataBase.area, INDataBase.dt, INDataBase.maxit,
                INDataBase.sources.STFLst, INDataBase.receivers.StaLst, INDataBase.snapshots, INDataBase.padding);
            # This part makes an animation using matplotlibs animation API
            background = (INDataBase.velocityM.velocity - INDataBase.velocityM.Vback)/ INDataBase.velocityM.Vback
            fig = mpl.figure(figsize=(16, 12))
            if len(INDataBase.PlotRec) != 0:
                Nplot=len(INDataBase.PlotRec)+1
                ax=mpl.subplot(Nplot,1,1)
                # ax = plt.subplot2grid((Nplot,1),(0, 2), rowspan=2)
                # ax=mpl.subplots()
                ax.xaxis.tick_top()
                ax.xaxis.set_label_position('top') 
                for i in INDataBase.PlotRec:
                    mpl.points(INDataBase.receivers.StaLst[i], '^b', size=8);
            # mpl.subplots_adjust(right=0.98, left=0.11, hspace=0.5, top=0.93)
            if not (background.max()<0.0001 and background.min() > -0.0001):
                Vbackground = plt.contour(background, extent=INDataBase.area,
                                       cmap=mpl.cm.seismic_r, vmin=-1, vmax=1)
            wavefield= mpl.imshow(np.zeros_like(INDataBase.velocityM.velocity), extent=INDataBase.area,
                                   cmap=mpl.cm.seismic_r, vmin=vmin, vmax=vmax);
            # mpl.points(INDataBase.receivers.StaLst, '^b', size=8);
            mpl.points(INDataBase.sources.EventLst, '*g', size=8)
            mpl.ylim(INDataBase.area[2:][::-1])
            mpl.xlabel('x (km)')
            mpl.ylabel('z (km)')
            mpl.m2km()
            if len(INDataBase.PlotRec) != 0:
                # ax=mpl.subplots()
                global Plot_seismograms
                global seis_index
                seis_index=INDataBase.PlotRec;
                Plot_seismograms=[];
                global times
                times = np.linspace(0, INDataBase.dt * INDataBase.maxit, INDataBase.maxit)
                for i in range(len(INDataBase.PlotRec)):
                    mpl.subplot(Nplot,1,i+2)
                    # plt.subplot2grid((Nplot,1),(0, 2), rowspan=2)
                    seismogram1, = mpl.plot([], [], '-k', lw=3);
                    Plot_seismograms.append(seismogram1)
                    mpl.xlim(0, INDataBase.duration)
                    mpl.ylim(vmin/2, vmax/2)
                    mpl.title('X: '+str(INDataBase.receivers.StaLst[INDataBase.PlotRec[i]][0]/1000.)+'km Z:'+
                              str(INDataBase.receivers.StaLst[INDataBase.PlotRec[i]][1]/1000.)+'km')
                ani = animation.FuncAnimation( fig, animate1,  frames=INDataBase.maxit / INDataBase.snapshots, interval=1)
            else:
                ani = animation.FuncAnimation( fig, animate2,  interval=1)
            if outmp4!='':
                ani.save(outmp4,dpi=300)
            mpl.show()
            # self.seismograms = seismograms;
        return
    
    def SaveSeismograms(self, outdir, PRX='MEM2D_'):
        # try:
        print 'ATTENTION: Saving seismograms!'
        INDataBase = self.INDataBase;
        dt=INDataBase.dt;
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        for i in np.arange(len(INDataBase.receivers.StaLst)):
            X=INDataBase.receivers.StaLst[i][0];
            Z=INDataBase.receivers.StaLst[i][1];
            tr=obspy.core.Trace();
            # Data
            tr.data=self.seismograms[i];
            # Header
            tr.stats['sac']={};
            tr.stats.sac.stlo=X/1000.;
            tr.stats.sac.stla=Z/1000.;
            tr.stats.sac.user0=INDataBase.receivers.ds;
            Nx=int(X/INDataBase.receivers.ds);
            Nz=int(Z/INDataBase.receivers.ds);
            tr.stats.station=str(Nx)+'_'+str(Nz);
            tr.stats.delta=dt;
            tr.stats.network='MEM2D';
            tr.stats.sac.b=0;
            tr.stats.sac.evlo=INDataBase.sources.EventLst[0][0]/1000.;
            tr.stats.sac.evla=INDataBase.sources.EventLst[0][1]/1000.;
            tr.stats.sac.dist=np.sqrt( (tr.stats.sac.stlo-tr.stats.sac.evlo)**2 + (tr.stats.sac.stla-tr.stats.sac.evla)**2 );
            sacfname=outdir +'/'+ PRX + tr.stats.station +'.SAC';
            tr.write(sacfname, format='sac');
        return;
            
    
     
def animate1(i):
    t, u, seismograms = simulation.next()
    for i in np.arange(len(Plot_seismograms)):
        Plot_seismograms[i].set_data(times[:t + 1], seismograms[seis_index[i]][:t + 1]) 
    wavefield.set_array(u[::-1])
    return 
        
    
def animate2(i):
    t, u, seismograms = simulation.next()
    # for i in np.arange(len(Plot_seismograms)):
    #     Plot_seismograms[i].set_data(times[:t + 1], seismogram[i][:t + 1]) 
    wavefield.set_array(u[::-1])
    return 
        

    
        

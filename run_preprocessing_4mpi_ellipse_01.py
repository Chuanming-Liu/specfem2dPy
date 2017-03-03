#!/usr/bin/env python
import stations
import numpy as np
import vmodel


### Generate Station List
# SLst=stations.StaLst()
# # # # SLst.HomoStaLst(xmin=0, xmax=16000000, dx=20000, zmin=500000, zmax=500000, dz=20000)
# SLst.HomoStaLst(xmin=0, xmax=8000000, dx=20000, zmin=2000000, zmax=2000000, dz=20000)
# SLst.WriteStaList('/lustre/janus_scratch/life9360/specfem2d_working_dir/cosine_anomaly_8000km_4000km_R_100km/DATA/STATIONS')

### Velocity Model

# # datadir = '/lustre/janus_scratch/life9360/specfem2d_working_dir/cosine_anomaly_16000km_1000km_R_400km/DATA/model'
# # outdir ='/lustre/janus_scratch/life9360/specfem2d_working_dir/cosine_anomaly_16000km_1000km_R_400km/DATA'
# 
# datadir = '/lustre/janus_scratch/life9360/specfem2d_working_dir/cosine_anomaly_8000km_4000km_R_350km/DATA/model'
# outdir ='/lustre/janus_scratch/life9360/specfem2d_working_dir/cosine_anomaly_8000km_4000km_R_350km/DATA'


datadir = '/lustre/janus_scratch/life9360/specfem2d_working_dir/homo_lens_vert_8000km_4000km_R_250km_d_100km/DATA/model'
outdir ='/lustre/janus_scratch/life9360/specfem2d_working_dir/homo_lens_vert_8000km_4000km_R_250km_d_100km/DATA'
# 
print outdir
R = 250000.
nproc = 12
for iproc in xrange(nproc):
    print iproc
    infname = datadir + '/proc%06d' % iproc + '_rho_vp_vs.dat'
    outfname = outdir + '/proc%06d' % iproc + '_rho_vp_vs.dat'
    Vm=vmodel.vmodel(xmin=0, xmax=8000000, Nx=1600, zmin=0, zmax=4000000, Nz=800, Vs=3000.)
    # Vm=vmodel.vmodel(xmin=0, xmax=10000000, Nx=2000, zmin=0, zmax=2000000, Nz=400, Vs=3000.)
    Vm.read(infname)
    Vm.setbackground(vs=3000.)
    # Vm.LensHomoModel(z0=1000000, x0=3500000, d1=50000, d2=50000, r1=R, r2=R, zmin=None, zmax=None, va=None, dv=-0.1)
    Vm.LensHomoModel_v2(z0=1000000, x0=3500000, d1=50000, d2=50000, r1=R, r2=R, zmin=None, zmax=None, va=None, dv=-0.1)
    # Vm.CircleCosineAnomaly( Xc=3500000, Zc=2000000, R=350000, dv = -0.1)
    # Vm.BlockHomoAnomaly( Xmin=1400000, Xmax=1600000, Zmin=0, Zmax=2000000, dv = -0.1)
    Vm.write(outfname, dt=0.05, fc=0.1)
    del Vm
    
# Vm=vmodel.vmodel(xmin=0, xmax=8000000, Nx=1600, zmin=0, zmax=4000000, Nz=800, Vs=3000.)
# Vm.LensHomoModel_v2(z0=2000000, x0=3500000, d1=50000, d2=50000, r1=300000, r2=300000, zmin=None, zmax=None, va=None, dv=-0.1)
# # # # Vm.CircleHomoAnomaly( Xc=1600000, Zc=1000000, R=100000, va = 2989.6)
# # # Vm.BlockHomoAnomaly( Xmin=1400000, Xmax=1600000, Zmin=0, Zmax=2000000, dv = -0.1)
# Vm.plot(vmin=2.7, vmax=3.3, cmap='seismic_r')

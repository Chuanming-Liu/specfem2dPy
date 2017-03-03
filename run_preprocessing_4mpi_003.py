#!/usr/bin/env python
import stations
import numpy as np
import vmodel

Dx=1000.*1000.
Dz=1000.*700.

### Generate Station List
# SLst=stations.StaLst()
# # SLst.HomoStaLst(xmin=0, xmax=3000000+2*Dx, dx=20000, zmin=300000+Dz, zmax=300000+Dz, dz=20000)
# SLst.HomoStaLst(xmin=0, xmax=6000000+2*Dx, dx=20000, zmin=300000+Dz, zmax=300000+Dz, dz=20000)
# SLst.WriteStaList('/lustre/janus_scratch/life9360/specfem2d_working_dir/2d_3d_staircase_basin/DATA/STATIONS')

### Velocity Model
# datadir = '/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing/DATA/model'
# outdir ='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing/DATA'
# # 
datadir = '/lustre/janus_scratch/life9360/specfem2d_working_dir/2d_3d_staircase_basin_8000km_2000km_R_100km/DATA/model'
outdir ='/lustre/janus_scratch/life9360/specfem2d_working_dir/2d_3d_staircase_basin_8000km_2000km_R_100km/DATA'
# 
# Vm=vmodel.vmodel(xmin=0, xmax=3000000+2*Dx, Nx=800, zmin=0, zmax=600000+2*Dz, Nz=400, Vs=3231.62)
# Vm.ASDFmodel('/lustre/janus_scratch/life9360/specfem2d_working_dir/2d_3d_ring_basin/eigendisp_ringbasin.h5', Dx=Dx, Dz=Dz)
# Vm.plot(vmin=2.7, vmax=3.3)
# 
# nproc = 12
# for iproc in xrange(nproc):
#     print iproc
#     infname = datadir + '/proc%06d' % iproc + '_rho_vp_vs.dat'
#     outfname = outdir + '/proc%06d' % iproc + '_rho_vp_vs.dat'
#     # Vm=vmodel.vmodel(xmin=0, xmax=3000000+2*Dx, Nx=800, zmin=0, zmax=600000+2*Dz, Nz=400, Vs=3231.62)
#     # Vm.read(infname)
#     # Vm.ASDFmodel('/lustre/janus_scratch/life9360/cps_working_dir/eigendisp_ringbasin.h5', Dx=Dx, Dz=Dz)
#     # Vm=vmodel.vmodel(xmin=0, xmax=4000000, Nx=800, zmin=0, zmax=2000000, Nz=400, Vs=3000.)
#     Vm=vmodel.vmodel(xmin=0, xmax=8000000, Nx=1600, zmin=0, zmax=2000000, Nz=400, Vs=3000.)
#     Vm.read(infname)
#     Vm.setbackground(vs=3231.4)
#     Vm.CircleHomoAnomaly( Xc=2100000, Zc=1000000, R=100000, va = 2989.6)
#     # Vm.CircleHomoAnomaly( Xc=1600000, Zc=1000000, R=100000, va = 2989.6)
#     Vm.write(outfname, dt=0.05, fc=0.1)
#     del Vm
    
Vm=vmodel.vmodel(xmin=0, xmax=4000000, Nx=800, zmin=0, zmax=2000000, Nz=400, Vs=3231.4)
Vm.CircleHomoAnomaly( Xc=1600000, Zc=1000000, R=100000, va = 2989.6)
Vm.plot(vmin=2.9896, vmax=3.2314, cmap='YlOrRd_r')

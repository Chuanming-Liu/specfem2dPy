#!/usr/bin/env python
import stations
import numpy as np
import vmodel

### Generate Station List
# SLst=stations.StaLst()
# SLst.HomoStaLst(xmin=0, xmax=8000000, dx=10000, zmin=2000000, zmax=2000000, dz=10000)
# SLst.WriteStaList('/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing/DATA/STATIONS')

### Velocity Model
datadir = '/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing/DATA/model'
outdir ='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing/DATA'
# 
# datadir = '/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_EA_001/DATA/model_EA'
# outdir ='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_EA_001/DATA'
# 
# nproc = 12
# for iproc in xrange(nproc):
#     print iproc
#     infname = datadir + '/proc%06d' % iproc + '_rho_vp_vs.dat'
#     outfname = outdir + '/proc%06d' % iproc + '_rho_vp_vs.dat'
#     # Vm=vmodel.vmodel(xmin=0, xmax=2000000, Nx=400, zmin=0, zmax=2000000, Nz=400, Vs=3000.)
#     Vm=vmodel.vmodel(xmin=0, xmax=8000000, Nx=1600, zmin=0, zmax=4000000, Nz=800, Vs=3000.)
#     Vm.read(infname)
#     Vm.setbackground(vs=3000.)
#     Vm.CircleCosineAnomaly( Xc=1500000, Zc=2000000, R=100000, dv = -0.1)
#     # Vm.readPhv(x0=0, z0=0, infname='10.phase_smooth.map')
#     Vm.write(outfname, dt=0.05, fc=0.1)
#     del Vm
    
Vm=vmodel.vmodel(xmin=0, xmax=8000000, Nx=1600, zmin=0, zmax=2000000, Nz=400, Vs=3000.)
Vm.CircleCosineAnomaly( Xc=1500000, Zc=1000000, R=100000, dv = -0.1)
Vm.plot(vmin=2.7, vmax=3.3)


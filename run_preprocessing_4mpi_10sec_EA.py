#!/usr/bin/env python
import stations
import numpy as np
import vmodel

### Generate Station List
# SLst=stations.StaLst()
# SLst.HomoStaLst(xmin=0, xmax=4000000, dx=20000, zmin=0, zmax=1000000, dz=20000)
# SLst.WriteStaList('/lustre/janus_scratch/life9360/specfem2d_working_dir/kernel_travelT/DATA/STATIONS')

### Velocity Model
datadir = '/lustre/janus_scratch/life9360/specfem2d_working_dir/kernel_travelT_0.003/DATA/model'
outdir ='/lustre/janus_scratch/life9360/specfem2d_working_dir/kernel_travelT_0.003/DATA'
# 
# datadir = '/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_EA_001/DATA/model_EA'
# outdir ='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_EA_001/DATA'
#


nproc = 12
for iproc in xrange(nproc):
    print iproc
    infname = datadir + '/proc%06d' % iproc + '_rho_vp_vs.dat'
    outfname = outdir + '/proc%06d' % iproc + '_rho_vp_vs.dat'
    # Vm=vmodel.vmodel(xmin=0, xmax=2000000, Nx=400, zmin=0, zmax=2000000, Nz=400, Vs=3000.)
    Vm=vmodel.vmodel(xmin=0, xmax=4000000, Nx=800, zmin=0, zmax=2000000, Nz=400, Vs=3000.)
    Vm.read(infname)
    Vm.setbackground(vs=3000.)
    Vm.CircleCosineAnomaly( Xc=1000000, Zc=1000000, R=100000, dv = -0.003)
    # Vm.BlockHomoAnomaly( Xmin=900000, Xmax=1100000, Zmin=0, Zmax=2000000, dv = -0.1)
    # Vm.readPhv(x0=0, z0=0, infname='10.phase_smooth.map')
    Vm.write(outfname, dt=0.05, fc=0.1)
    del Vm

### plot figure
# Vm=vmodel.vmodel(xmin=0, xmax=4000000, Nx=800, zmin=0, zmax=2000000, Nz=400, Vs=3000.)
# Vm.CircleCosineAnomaly( Xc=1000000, Zc=1000000, R=100000, dv = -0.1)
# Vm.plot(vmin=2.7, vmax=3.3)
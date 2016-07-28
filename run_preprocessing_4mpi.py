#!/usr/bin/env python
import stations
import numpy as np
import vmodel

### Generate Station List
# SLst=stations.StaLst()
# SLst.HomoStaLst(xmin=0, xmax=2000000, dx=20000, zmin=0, zmax=2000000, dz=20000)
# SLst.WriteStaList('STATIONS')

### Velocity Model
datadir = '/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_D_570/DATA/model'
outdir ='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_D_570/DATA'
nproc = 12
for iproc in xrange(nproc):
    infname = datadir + '/proc%06d' % iproc + '_rho_vp_vs.dat'
    outfname = outdir + '/proc%06d' % iproc + '_rho_vp_vs.dat'
    Vm=vmodel.vmodel(xmin=0, xmax=2000000, Nx=400, zmin=0, zmax=2000000, Nz=400, Vs=3000.)
    Vm.read(infname)
    Vm.setbackground(vs=3000.)
    Vm.CircleCosineAnomaly( Xc=1200000, Zc=1200000, R=100000, dv = -0.1)
    Vm.write(outfname, dt=0.05, fc=0.1)
    del Vm


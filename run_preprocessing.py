#!/usr/bin/env python
import stations
import numpy as np
import vmodel


### Generate Station List
# SLst=stations.StaLst()
# SLst.HomoStaLst(xmin=0, xmax=2000000, dx=20000, zmin=0, zmax=2000000, dz=20000)
# SLst.WriteStaList('STATIONS')

### Velocity Model
# Vm=vmodel.vmodel(xmin=0, xmax=2000000, Nx=400, zmin=0, zmax=2000000, Nz=400, Vs=3000.)
Vm=vmodel.vmodel(xmin=0, xmax=8000000, Nx=1600, zmin=0, zmax=4000000, Nz=800, Vs=3000.)
# Vm.read('/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_D/DATA/proc000000_rho_vp_vs.dat_backup')
# Vm.read('/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_D_1200/DATA/proc000000_rho_vp_vs.dat_backup')
Vm.setbackground(vs=3000.)
Vm.readPhv(x0=200000, z0=200000, evlo=129., evla=41.306, infname='10.phase_smooth.map')
# Vm.CircleCosineAnomaly( Xc=1200000, Zc=1200000, R=100000, dv = -0.1)
Vm.plot()
# Vm.write('/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_D/DATA/proc000000_rho_vp_vs.dat', dt=0.05, fc=0.1)
# Vm.write('./test.dat', dt=0.05, fc=0.1)
### Velocity Model
# Vm=vmodel.vmodel(xmin=0, xmax=3000000+2*Dx, Nx=640, zmin=0, zmax=600000+2*Dz, Nz=520, Vs=3023.22)
# Vm.read('/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_0.1_20/DATA/proc000000_rho_vp_vs.dat')
# Vm.CircleHomoAnomaly(Xc=1500000+Dx+800000, Zc=300000+Dz, R=100000, va=3395.37)
# Vm.CircleHomoAnomaly(Xc=1500000+Dx-800000, Zc=300000+Dz, R=100000, va=2683.98)
# Vm.write('/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_0.1_20/DATA/proc000000_rho_vp_vs.dat')
# # 
# Vm.write('/projects/life9360/code/SEM/specfem2d/EXAMPLES/LFMembrane_SH_0.1_20/DATA/proc000000_rho_vp_vs.dat',
#          dt=0.05, fc=0.1)
# Vm.plot()
# Vm.write('proc000000_rho_vp_vs.dat', dt=0.05, fc=0.1)

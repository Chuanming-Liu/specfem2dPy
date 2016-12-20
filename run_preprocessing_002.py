#!/usr/bin/env python
import stations
import numpy as np
import vmodel

Dx=1000.*100.
Dz=1000.*1000.

### Generate Station List
SLst=stations.StaLst();
SLst.HomoStaLst(xmin=0, xmax=1000000, zmin=0, zmax=1000000, dx=10000, dz=10000);
SLst.WriteStaList('/projects/life9360/code/SEM/specfem2d/EXAMPLES/LFMembrane_SH/DATA/STATIONS');
# 
# ### Velocity Model
Vm=vmodel.vmodel(xmin=0, xmax=3000000+2*Dx, Nx=640, zmin=0, zmax=600000+2*Dz, Nz=520, Vs=3231.62)
Vm.ASDFmodel('../CPSPy/eigendisp_ringbasin.h5')
Vm.plot()
# # Vm.CircleGradualAnomaly(Xc=640000, Zc=640000, R=60000, Va=2500);
# 
# Vm.read('/projects/life9360/code/SEM/specfem2d/EXAMPLES/LFMembrane_SH/DATA/proc000000_rho_vp_vs.dat')
# Vm.CircleCosineAnomaly(Xc=700000, Zc=700000, R=60000, Va=2500);
# # # # # Vm.BlockHomoAnomaly(Xmin=50000, Xmax=100000, Zmin=120000, Zmax=180000, Va=3000);
# # # # # Vm.CircleHomoAnomaly(Xc=100000, Zc=100000, R=50000, Va=2000.)
# # # # # Vm.BlockHomoAnomaly(Xmin=0, Xmax=100000, Zmin=0, Zmax=180000, Va=3000);
# Vm.plot()
# Vm.write('/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH/DATA/model.dat', dt=0.05, fc=0.1)
# # # IC=specfem2dPy.InputChecker(dt=0.01, dx=5., dz=5., fc=0.1, lpd=4, vmin=3.0, vmax=3.5)
# 
# itest=np.loadtxt('/projects/life9360/code/SEM/specfem2d/EXAMPLES/LFMembrane/DATA/proc000000_rho_vp_vs.dat');
# AAA=itest[:,0];
# BBB=itest[:,1];

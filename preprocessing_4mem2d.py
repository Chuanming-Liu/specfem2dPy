#!/usr/bin/env python
import stations
import numpy as np
import vmodel


### Generate Station List
SLst=stations.StaLst();
SLst.HomoStaLst(xmin=20000, xmax=2980000, zmin=20000, zmax=580000, dx=20000, dz=20000)
SLst.WriteStaList('/projects/life9360/code/SEM/specfem2d/EXAMPLES/LFMembrane_SH_0.1_20/DATA/STATIONS');
# 
### Velocity Model
Vm=vmodel.vmodel(xmin=0, xmax=3000000, Nx=600, zmin=0, zmax=600000, Nz=120, Vs=3231.62)
# Vm.CircleGradualAnomaly(Xc=640000, Zc=640000, R=60000, Va=2500);
# 
Vm.write('/projects/life9360/code/SEM/specfem2d/EXAMPLES/LFMembrane_SH_0.1_20/DATA/proc000000_rho_vp_vs.dat',
         dt=0.05, fc=0.1)


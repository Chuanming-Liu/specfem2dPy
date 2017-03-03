#!/usr/bin/env python
import stations
import numpy as np
import vmodel

Dx=1000.*1000.
Dz=1000.*700.


    
Vm=vmodel.vmodel(xmin=0, xmax=4000000, Nx=100, zmin=0, zmax=2000000, Nz=50, Vs=3231.4)
Vm.CircleHomoAnomaly( Xc=1600000, Zc=1000000, R=100000, va = 2989.6)
Vm.plot(vmin=2.9896, vmax=3.2314, cmap='YlOrRd_r')

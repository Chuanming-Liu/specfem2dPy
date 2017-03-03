#!/usr/bin/env python
import specfem2dPy
import numpy as np
# import symData2d
import field2d_cartesian as field2d
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy.ndimage.filters import convolve, gaussian_filter

# Mask=np.ones([15,15])/225;
# def Convolve(im, inMask):
#     return convolve(ma.getdata(im), inMask);


### Field Analysis
dx=20.
dy=20.
# Tfield=field2d.Field2d(Nx=1000/dx, Ny=1000/dx, dx=dx, dy=dy);
# dx=5.
# dy=5.
Tfield=field2d.Field2d(Nx=2000/dx, Ny=2000/dx, dx=dx, dy=dy);
# Tfield.LoadFile('10km/Tph_10.0.txt')
Tfield.LoadFile('./D_990km/Amp_10.0.txt')
Tfield.natgridInterp();
Tfield.PlotField()
# Tfield.CuttingEdges(nx=20, ny=20, fieldtype='TravelT')
# Tfield.Gradient()
# Tfield.GetDeflection()
# vabs=10
# Tfield.PlotDeflection(vmin=-vabs, vmax=vabs)

# Tfield.get_polar_value_deflection(mindist=100, ndist=7, ddist=100., dper=10, refangle=45, plotflag=True)
# Tfield.get_polar_value_deflection(distLst=[700., 1000., 1200., 1450., 1800.], dper=10, refangle=45, minpolar=-15, maxpolar=15, plotflag=True)
# Tfield.get_polar_value_deflection(distLst=[300, 700., 900., 1300., 1700.], dper=10, refangle=45, minpolar=-20, maxpolar=20, plotflag=True)
# Tfield.get_polar_value_deflection(distLst=[200, 400., 700., 1150., 1500.], dper=10, refangle=45, minpolar=-25, maxpolar=25, plotflag=True)
# Tfield.get_polar_value_deflection(distLst=[150., 300., 500., 800.], dper=10, refangle=45, minpolar=-80, maxpolar=75, plotflag=True)
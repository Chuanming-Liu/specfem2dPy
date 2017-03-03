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
Afield=field2d.Field2d(Nx=2000/dx, Ny=2000/dx, dx=dx, dy=dy);
# Tfield.LoadFile('10km/Tph_10.0.txt')
Afield.LoadFile('./D_1200km/Amp_10.0.txt')
Afield.natgridInterp();

# Tfield.CuttingEdges(nx=20, ny=20, fieldtype='TravelT')
# Tfield.Gradient()
# Tfield.GetDeflection()
# vabs=10
# Tfield.PlotDeflection(vmin=-vabs, vmax=vabs)

Afield.get_polar_value(distLst=[1100., 1200., 1300., 1400., 1500.], dper=8, refangle=45, minpolar=-30, maxpolar=30, plotflag=True)
Afield.PlotField(contourflag=False)



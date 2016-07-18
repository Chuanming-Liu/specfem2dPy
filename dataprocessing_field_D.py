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
Tfield.LoadFile('./D_670km/Tph_10.0.txt')
Tfield.natgridInterp();
# Tfield.PlotField()
# Tfield.CuttingEdges(nx=20, ny=20, fieldtype='TravelT')
Tfield.Gradient()
Tfield.GetDeflection()
vabs=6
Tfield.PlotDeflection(vmin=-vabs, vmax=vabs)


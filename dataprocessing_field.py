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

# 
# ### Field Analysis
# dx=10.
# dy=10.
# # Tfield=field2d.Field2d(Nx=1000/dx, Ny=1000/dx, dx=dx, dy=dy);
# # dx=5.
# # dy=5.
# Tfield=field2d.Field2d(Nx=480/dx, Ny=480/dx, dx=dx, dy=dy);
# # Tfield.LoadFile('10km/Tph_10.0.txt')
# Tfield.LoadFile('./field_homo/TravelT.gr.10.0.txt')
# Tfield.XarrIn = Tfield.XarrIn /1000.
# Tfield.YarrIn = Tfield.YarrIn / 1000.
# Tfield.natgridInterp();
# # Tfield.PlotField()
# # Tfield.CuttingEdges(nx=20, ny=20, fieldtype='TravelT')
# Tfield.Gradient()
# Tfield.GetDeflection()
# Tfield.PlotDeflection(vmin=-5, vmax=5)
# # Tfield.PlotPropagation()
# # Tfield.Gradient(method='convolve', order=2);
# # Tfield.GetApparentV();
# # Tfield.PlotAppV()
# # Tfield.GetDoT();
# # Tfield.PlotDoT(vmin = 3.4, vmax = 3.6);
# # plt.show()
# ### Travel Time Laplacian
# # # # Tfield.LaplacianEqualXY();
# # # # Tfield.PlotLaplacian();
# # # # Tfield.Laplacian();
# # # # Tfield.PlotLaplacian();
# # # # Tfield.Laplacian(method='convolve',order=2);
# # # # Tfield.PlotLaplacian();
# # # # Tfield.Laplacian(method='convolve',order=4);
# # # # Tfield.PlotLaplacian();
# # # # Tfield.Laplacian(method='convolve',order=6);
# # # # Tfield.PlotLaplacian();
# 
# 
dx=10.
dy=10.
Afield=field2d.Field2d(Nx=1000/dx, Ny=1000/dx, dx=dx, dy=dy);
# dx=5.
# dy=5.
# Tfield=field2d.Field2d(Nx=480/dx, Ny=480/dx, dx=dx, dy=dy);
Afield.LoadFile('field_circle/Amplitude.10.0.txt');
Afield.natgridInterp();
Afield.Laplacian(method='convolve',order=2);
Afield.GetLplcCorrection(10.);
Afield.PlotLplcCo();
Tfield.GetCorV(Afield);
Tfield.PlotCorV();
plt.show()

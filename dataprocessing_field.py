#!/usr/bin/env python
import specfem2dPy
import numpy as np
import symData2d
import field2d_cartesian as field2d
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy.ndimage.filters import convolve, gaussian_filter

# Mask=np.ones([15,15])/225;
# def Convolve(im, inMask):
#     return convolve(ma.getdata(im), inMask);


### Field Analysis
dx=5.
dy=5.
myfield=field2d.Field2d(Nx=1000/dx, Ny=1000/dx, dx=dx, dy=dy);
myfield.LoadFile('/lustre/janus_scratch/life9360/specfem2d_data/field_homo/TravelT.ph.10.0.txt')
myfield.natgridInterp();
# myfield.PlotField()
# myfield.CuttingEdges(nx=20, ny=20, fieldtype='TravelT')
# myfield.Gradient();
myfield.Gradient(method='convolve', order=4);
myfield.GetApparentV();
myfield.PlotAppV()
myfield.GetDoT(enx=250., eny=250.);
myfield.PlotDoT();
# 
# dx=3.
# dy=3.
# myfield2=field2d.Field2d(Nx=480/dx, Ny=480/dx, dx=dx, dy=dy);
# myfield2.LoadFile('/lustre/janus_scratch/life9360/specfem2d_data/field_homo/TravelT.ph.10.0.txt')
# myfield2.natgridInterp();
# myfield2.CuttingEdges(nx=10, ny=10, fieldtype='TravelT')
# # myfield2.PlotField()
# myfield2.Gradient(method='convolve', order=6);
# myfield2.GetApparentV();
# myfield2.PlotAppV()
# 
# myfield2.Gradient(method='freq');
# # myfield2.Gradient(method='convolve',order=4);
# myfield2.GetApparentV();
# 
# # myfield2.AppV=Convolve(myfield2.AppV, Mask)
# # myfield2.AppV=gaussian_filter(myfield2.AppV, 4.)
# myfield2.PlotAppV()
# plt.show()
# # myfield2.Gradient()
# # grad=myfield2.grad[1]
# # outdiff1=myfield2.fftDiff(1,0)
# # outdiff2=myfield2.fftDiff2(1,0)
# # outdiff2=myfield2.fftDiff2(1,0)      
# # outdiff=myfield2.fftDiff2(1,1)
# 
# 
# # plt.subplot(131),plt.imshow(outdiff1, cmap = 'seismic',vmin=-0.1,vmax=0.1)
# # # plt.title('Input Image'), plt.xticks([]), plt.yticks([])
# # plt.subplot(132),plt.imshow(outdiff2, cmap = 'seismic',vmin=-0.1,vmax=0.1)
# # plt.subplot(133),plt.imshow(outdiff1-outdiff2, cmap = 'seismic',vmin=-0.1,vmax=0.1)
# # plt.colorbar()
# # # plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])
plt.show()

# plt.show()
# myfield.LaplacianEqualXY()
# myfield.GetLplcCorrection(per=4.0)
# # myfield.Laplacian()
# # myfield.PlotLplcCo()

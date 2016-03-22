#!/usr/bin/env python
import specfem2dPy
import numpy as np
import symData2d
import field2d_cartesian as field2d
import numpy.ma as ma
import matplotlib.pyplot as plt
from skimage.filters import roberts, sobel_v, scharr_h, scharr_v, prewitt
def makeGaussian(size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """
    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    X,Y=np.meshgrid(x,y)
    return X,Y, np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)


dx=3.
dy=3.
myfield2=field2d.Field2d(Nx=480/dx, Ny=480/dx, dx=dx, dy=dy);
X,Y, myfield2.Zarr=makeGaussian(161, fwhm=50);
outdiff1=myfield2.fftDiff(1,0)
outdiff2=myfield2.fftDiff2(1,0,1024,1024)      
edge_sobel = sobel_v(myfield2.Zarr)


plt.subplot(131),plt.imshow(outdiff1, cmap = 'seismic',vmin=-0.01,vmax=0.01)
# plt.title('Input Image'), plt.xticks([]), plt.yticks([])
plt.subplot(132),plt.imshow(outdiff2, cmap = 'seismic',vmin=-0.01,vmax=0.01)
plt.subplot(133),plt.imshow(edge_sobel/5, cmap = 'seismic',vmin=-0.01,vmax=0.01)
plt.colorbar()
# plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])
plt.show()


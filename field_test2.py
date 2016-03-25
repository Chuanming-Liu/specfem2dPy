#!/usr/bin/env python
import specfem2dPy
import numpy as np
import symData2d
import field2d_cartesian as field2d
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy.ndimage import convolve, binary_erosion, generate_binary_structure
X_diff_weight_2 = np.array([1., 0., -1.])/2.;
Y_diff_weight_2 = X_diff_weight_2.T;

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
myfield=field2d.Field2d(Nx=480/dx, Ny=480/dx, dx=dx, dy=dy);
X,Y, myfield.Zarr=makeGaussian(161, fwhm=100);
myfield.Gradient()
outdiff1=myfield.grad[0];
myfield.Gradient(method='convolve', order=2);
outdiff2=myfield.grad[0];
myfield.Gradient(method='convolve', order=4);
outdiff3=myfield.grad[0];
myfield.Gradient(method='convolve', order=6);
outdiff4=myfield.grad[0];
myfield.Gradient(method='freq');
outdiff5=myfield.grad[0];

# plt.subplot(141),plt.imshow(outdiff1, cmap = 'seismic',vmin=-0.01,vmax=0.01)
# plt.subplot(142),plt.imshow(outdiff2, cmap = 'seismic',vmin=-0.01,vmax=0.01)
# plt.subplot(143),plt.imshow(outdiff3, cmap = 'seismic',vmin=-0.01,vmax=0.01)
# plt.subplot(144),plt.imshow(outdiff5, cmap = 'seismic',vmin=-0.01,vmax=0.01)

# plt.subplot(141),plt.imshow(outdiff1, cmap = 'seismic',vmin=-0.01,vmax=0.01)
# plt.subplot(142),plt.imshow(outdiff2, cmap = 'seismic',vmin=-0.01,vmax=0.01)
# plt.subplot(143),plt.imshow(outdiff3, cmap = 'seismic',vmin=-0.01,vmax=0.01)
# plt.subplot(144),plt.imshow(outdiff5, cmap = 'seismic',vmin=-0.01,vmax=0.01)
plt.figure()
plt.subplot(141),plt.imshow(outdiff5, cmap = 'seismic',vmin=-0.01,vmax=0.01)
plt.colorbar()
plt.subplot(142),plt.imshow(outdiff2-outdiff5, cmap = 'seismic',vmin=-0.00001,vmax=0.00001)
plt.colorbar()
plt.subplot(143),plt.imshow(outdiff3-outdiff5, cmap = 'seismic',vmin=-0.00001,vmax=0.00001)
plt.colorbar()
plt.subplot(144),plt.imshow(outdiff4-outdiff5, cmap = 'seismic',vmin=-0.00001,vmax=0.00001)
plt.colorbar()

# plt.colorbar()
# plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])
plt.show()


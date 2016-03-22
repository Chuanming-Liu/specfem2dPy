import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.ndimage.filters import gaussian_filter
from math import pi
# import cv2





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

X,Y, Z=makeGaussian(100, fwhm=10);
# img = cv2.imread('messi5.jpg',0)
h = np.fft.fft2(Z)
hf = np.fft.fftshift(h)
# magnitude_spectrum = np.abs(fshift)
# u=np.arange(100.);
# v=np.arange(100.);
u=np.arange(100.)-50.;
v=np.arange(100.)-50.;
U,V,=np.meshgrid(u,v);
m = 1; n = 1;
hf2 = ((1j*2*pi*U/100)**m)*((1j*2*pi*V/100)**n)*hf;
out_diff1 = np.real(np.fft.ifft2(np.fft.ifftshift(hf2)));

u=np.arange(100.)-49.;
v=np.arange(100.)-49.;
U,V,=np.meshgrid(u,v);
m = 1; n = 1;
hf2 = ((1j*2*pi*U/100)**m)*((1j*2*pi*V/100)**n)*hf;
out_diff2 = np.real(np.fft.ifft2(np.fft.ifftshift(hf2)));

# plt.subplot(121),plt.imshow(Z, cmap = 'seismic')
# plt.title('Input Image'), plt.xticks([]), plt.yticks([])
# plt.subplot(122),plt.imshow(out_diff, cmap = 'seismic')
# plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])
# plt.show()


# X,Y, Z=makeGaussian(100, fwhm=50);
# spectrum=np.fft.fft2(Z)
# # Z=gaussian_filter(np.ones([100,100]), sigma=100)
# 
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# # surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
# #                        linewidth=0, antialiased=False)
# surf = ax.plot_surface(X, Y, out_diff, rstride=1, cstride=1, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
# # 
# # ax.imshow(Z)
# plt.show()
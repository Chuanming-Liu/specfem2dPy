#!/usr/bin/env python
import numpy as np
import scipy.interpolate as sinter
import matplotlib.pyplot as plt

fname='/lustre/janus_scratch/yeti4009/for_lili/noisy.txt';
In=np.loadtxt(fname)
freq=In[:,0];
yp=In[:,1];
# yy=np.arange(1001.)/1000.;

# plt.plot(freq,y)
# plt.show()

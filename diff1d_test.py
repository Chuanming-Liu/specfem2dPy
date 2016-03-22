import numpy as np
import matplotlib.pyplot as plt

dt=0.5;
t=np.arange(100)*dt;
x=np.sin(t)
diffx=np.cos(t)
hf0=np.fft.fft(x);
hf1=np.fft.fft(x, n=128);
hf0=np.fft.fftshift(hf0);
hf1=np.fft.fftshift(hf1);
f0=np.arange(100)-50;
f1=np.arange(128)-64;
diff_hf0=1j*2.*np.pi*f0/100*hf0;
diff_hf1=1j*2.*np.pi*f1/128*hf1;

out_diff0=np.real( np.fft.ifft( np.fft.ifftshift(diff_hf0) ) )/dt;
out_diff1=np.real( np.fft.ifft( np.fft.ifftshift(diff_hf1) ) )/dt;
plt.plot(out_diff0,'x--')
plt.plot(out_diff1,'o-')
plt.plot(diffx,'^-')
import numpy as np
import matplotlib.pyplot as plt

dt=0.5;
t=np.arange(100)*dt;
x=np.sin(t)
hf0=np.fft.fft(x, n=64);
hf1=np.fft.fft(x, n=128);
hf2=np.fft.fft(x);
hf3=np.fft.fft(x,n=1024);
f0=1/dt*2*np.arange(64)/64
f1=1/dt*2*np.arange(128)/128
f2=1/dt*2*np.arange(100)/100
f3=1/dt*2*np.arange(1024)/1024
plt.figure()
plt.plot(f0, abs(hf0),'v-');
plt.plot(f1, abs(hf1),'x-');
plt.plot(f2, abs(hf2),'o--');
plt.plot(f3, abs(hf3),'^--');

ihf0=np.fft.ifft(hf0)
ihf1=np.fft.ifft(hf1)
ihf2=np.fft.ifft(hf2)
ihf3=np.fft.ifft(hf3)

plt.figure()
plt.plot(t,x,'k')
plt.plot(t[:64], ihf0,'v-');
plt.plot(np.arange(128)*dt, ihf1,'x--');
plt.plot(t, ihf2[:100],'o--');
plt.plot(np.arange(256)*dt, ihf3[:256],'^--');

plt.show()

# freq=1/dt*

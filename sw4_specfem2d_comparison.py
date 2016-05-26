import obspy
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


Dx=100.
# Dz=1000.*1000.
infname1 = './ak135_4mem2d_20km_0.1/Amp_10.0.txt'
infname2 = '../SW4Py/ak135_4mem2d_20km_0.1/Amp_10.0.txt'
inArr1=np.loadtxt(infname1)
inArr2=np.loadtxt(infname2)


XArr1=inArr1[:,0]-Dx
AmpArr1=inArr1[:,2]
ind1=np.argsort(XArr1)
XArr1=XArr1[ind1]
AmpArr1=AmpArr1[ind1]

XArr2=inArr2[:,0]
AmpArr2=inArr2[:,2]
ind2=np.argsort(XArr2)
XArr2=XArr2[ind2]
AmpArr2=AmpArr2[ind2]

XArr3= XArr2
AmpArr3= np.ones(XArr3.size) / np.sqrt(np.abs(XArr2-1500))

index1=np.where(XArr1==1300.)
amp1=AmpArr1[index1]
index2=np.where(XArr2==1300.)
amp2=AmpArr2[index2]
index3=np.where(XArr3==1300.)
amp3=AmpArr3[index3]

fig, ax=plt.subplots()
ax.plot(XArr1, AmpArr1/amp1,'o' );
ax.plot(XArr2, AmpArr2/amp2,'-x' );
ax.plot(XArr3, AmpArr3/amp3,'-' );
ax.fill_between(np.array([]), y1, y2, where=y2 >= y1, facecolor='green', interpolate=True)
plt.ylabel('Amplitude');
plt.xlabel('Distance(km)');

plt.show()

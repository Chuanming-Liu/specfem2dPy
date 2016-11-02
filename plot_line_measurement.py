import obspy
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


Dx=100.
# Dz=1000.*1000.
infname1 = '/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing/field_data/Amp_10.0.txt'
infname2 = '/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing/field_data/Amp_10.0.txt'
inArr1=np.loadtxt(infname1)
inArr2=np.loadtxt(infname2)


XArr1=inArr1[:,0]-Dx
AmpArr1=inArr1[:,2]
ind1=np.argsort(XArr1)
XArr1=XArr1[ind1]
AmpArr1=AmpArr1[ind1]
AmpArr1=AmpArr1[(XArr1<=1300) + (XArr1>=1700)]
XArr1=XArr1[(XArr1<=1300) + (XArr1>=1700)]

XArr2=inArr2[:,0]
AmpArr2=inArr2[:,2]
ind2=np.argsort(XArr2)
XArr2=XArr2[ind2]
AmpArr2=AmpArr2[ind2]
AmpArr2=AmpArr2[(XArr2<=1300) + (XArr2>=1700)]
XArr2=XArr2[(XArr2<=1300) + (XArr2>=1700)]

XArr3= XArr2
AmpArr3= np.ones(XArr3.size) / np.sqrt(np.abs(XArr2-1500))

index1=np.where(XArr1==1300.)
amp1=AmpArr1[index1]
index2=np.where(XArr2==1300.)
amp2=AmpArr2[index2]
index3=np.where(XArr3==1300.)
amp3=AmpArr3[index3]

fig, ax=plt.subplots()
ax.plot(XArr1, AmpArr1/amp1,'g-.o', lw=3, markersize=10, label='SPECFEM2D' );
ax.plot(XArr2, AmpArr2/amp2,'b-x', lw=3, markersize=10, label='SW4' );
ax.plot(XArr3, AmpArr3/amp3,'r-', lw=3, markersize=10, label='Predicted' );
ax.plot(1500., 1, 'y*', markersize=50)
plt.legend(loc='upper right', fontsize=25)
y1=600
y2=800
ax.fill_betweenx(np.array([0.2, 1.1]), y1, y2, facecolor='red', alpha=0.5)
y1=2200
y2=2400
ax.fill_betweenx(np.array([0.2, 1.1]), y1, y2, facecolor='blue', alpha=0.5)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
plt.ylabel('Amplitude', fontsize=30);
plt.xlabel('X position (km)', fontsize=30);
plt.title('Amplitude ', fontsize=30);
plt.show()

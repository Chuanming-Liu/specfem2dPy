import obspy
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats



Dx=500.
infname = '/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing/field_data/Tph_10.0.txt'
inArr=np.loadtxt(infname)

XArr=inArr[:,0]
TArr=inArr[:,2]

ind=np.argsort(XArr)
XArr=XArr[ind]
TArr=TArr[ind]

XArr2= XArr
TArr2= np.abs(XArr-Dx)/np.ones(XArr.size)/3.

# index=np.where(XArr==1300.)
# amp=AmpArr[index]
# index2=np.where(XArr2==1300.)
# amp2=AmpArr2[index2]
# index3=np.where(XArr3==1300.)
# amp3=AmpArr3[index3]

fig, ax=plt.subplots()
# ax.plot(XArr, TArr,'g-.o', lw=3, markersize=5, label='low' );
# ax.plot(XArr2, TArr2,'b--', lw=3, markersize=10, label='homo' );
ax.plot(XArr-Dx, (TArr-TArr2-1.4),'go', lw=3, markersize=10, label='travel time difference')
ax.plot(Dx-Dx, 1, 'y*', markersize=20)
# plt.legend(loc='upper right', fontsize=25, numpoints = 1)
y1=1400.-Dx
y2=1600-Dx
ax.fill_betweenx(np.array([-0.1, TArr.max()]), y1, y2, facecolor='red', alpha=0.5)
# y1=2200
# y2=2400
# ax.fill_betweenx(np.array([0.2, 1.1]), y1, y2, facecolor='blue', alpha=0.5)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
# plt.ylim([TArr.min(), TArr.max()])
plt.ylim([-0.1, (TArr-TArr2-1.4).max()])
plt.xlim([500., 7500.])
plt.xticks(np.arange(500., 8000., 500.))
plt.ylabel('Travel time difference (s)', fontsize=30);
plt.xlabel('Distance (km)', fontsize=30);
# plt.title('Travel Time with rectangle anomaly', fontsize=30);
plt.show()

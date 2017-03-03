import obspy
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def get_data(infname):
    Dx=500.
    inArr=np.loadtxt(infname)
    XArr=inArr[:,0]
    distArr=inArr[:,3]
    AmpArr=inArr[:,2]
    ind=np.argsort(XArr)
    XArr=XArr[ind]
    AmpArr=AmpArr[ind]
    distArr=distArr[ind]
    XArr2= XArr
    AmpArr2= np.ones(XArr.size) / np.sqrt(distArr)
    index=np.where(XArr==1000)
    amp=AmpArr[index]
    index2=np.where(XArr2==1000.)
    amp2=AmpArr2[index2]
    return XArr, AmpArr, amp, AmpArr2, amp2

fig, ax=plt.subplots()
Dx=500.
XArr, AmpArr, amp, AmpArr2, amp2=get_data(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/rectangle_4000km_2000km_D_200km/field_data/Amp_10.0.txt')
ax.plot(XArr-Dx, AmpArr/amp,'ro-', lw=3, markersize=15, label='rectangular anomaly' )
# XArr, AmpArr, amp, AmpArr2, amp2=get_data(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/cosine_anomaly_8000km_2000km_R_200km/field_data/Amp_10.0.txt')
# ax.plot(XArr-Dx, AmpArr/amp*amp2/AmpArr2,'b-', lw=3, markersize=10, label='R = 200 km' )
# XArr, AmpArr, amp, AmpArr2, amp2=get_data(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/cosine_anomaly_8000km_2000km_R_300km/field_data/Amp_10.0.txt')
# ax.plot(XArr-Dx, AmpArr/amp*amp2/AmpArr2,'r-', lw=3, markersize=10, label='R = 300 km' )
# XArr, AmpArr, amp, AmpArr2, amp2=get_data(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/cosine_anomaly_8000km_2000km_R_400km/field_data/Amp_10.0.txt')
# ax.plot(XArr-Dx, AmpArr/amp*amp2/AmpArr2,'g-', lw=3, markersize=10, label='R = 400 km' )
# XArr, AmpArr, amp, AmpArr2, amp2=get_data(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/cosine_anomaly_8000km_2000km_R_500km/field_data/Amp_10.0.txt')
# ax.plot(XArr-Dx, AmpArr/amp*amp2/AmpArr2,'m-', lw=3, markersize=10, label='R = 500 km' )

ax.plot(XArr-Dx, AmpArr2/amp2,'k--', lw=5, markersize=10, label='reference' );
plt.legend(loc='upper right', fontsize=20, numpoints = 1)
# y1=900.
# y2=1100.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='red', alpha=0.1)
# y1=800.
# y2=1200.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='red', alpha=0.1)
# y1=700.
# y2=1300.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='red', alpha=0.1)
# y1=600.
# y2=1400.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='red', alpha=0.1)
# y1=500.
# y2=1500.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='red', alpha=0.1)

y1=900.
y2=1100.
ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='r', alpha=0.5)
# y1=800.
# y2=900.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='b', alpha=0.5)
# y1=1100.
# y2=1200.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='b', alpha=0.5)
# y1=700.
# y2=800.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='r', alpha=0.5)
# y1=1200.
# y2=1300.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='r', alpha=0.5)
# y1=600.
# y2=700.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='g', alpha=0.5)
# y1=1300.
# y2=1400.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='g', alpha=0.5)
# y1=500.
# y2=600.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='m', alpha=0.5)
# y1=1400.
# y2=1500.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='m', alpha=0.5)

ax.tick_params(axis='x', labelsize=25)
ax.tick_params(axis='y', labelsize=25)
plt.yticks(np.arange(0.1, 1.5, 0.1))
plt.ylim([0.3, 1.5])
# plt.xticks(np.arange(200., 7700., 500.))
plt.xticks(np.arange(200., 7700., 500.))
plt.xlim([200., 3200.])
# plt.ylim([(TArr-TArr2-1.5).min(), (TArr-TArr2).max()])
# plt.ylabel('Normalized amplitude', fontsize=30);
plt.ylabel('Amplitude', fontsize=40);
# plt.xlabel('X position (km)', fontsize=30);
plt.xlabel('Distance (km)', fontsize=40);
# plt.title('Amplitude measurement with rectangle anomaly', fontsize=30);
# plt.title('Amplitude measurement with circular anomaly', fontsize=30);
plt.show()

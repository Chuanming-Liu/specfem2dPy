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

Aratio=np.array([]); flength=np.array([]); Amax=np.array([]); flength2=np.array([])
fig, ax=plt.subplots()
Dx=500.
XArr, AmpArr, amp, AmpArr2, amp2=get_data(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/homo_lens_10000km_2000km_R_100km_d_100km/field_data/Amp_10.0.txt')
ax.plot(XArr-Dx, AmpArr/amp,'k-', lw=5, markersize=10, label='R = 100 km' )
flength = np.append(flength, (XArr-Dx)[((AmpArr/amp)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp) ] - 3000.)
Amax = np.append(Amax, ((AmpArr/amp)[(XArr-Dx)> 3000.]).max())
flength2 = np.append(flength2, (XArr-Dx)[((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp*amp2/AmpArr2) ] - 3000.)
Aratio  = np.append(Aratio, ((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() )

XArr, AmpArr, amp, AmpArr2, amp2=get_data(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/homo_lens_10000km_2000km_R_150km_d_100km/field_data/Amp_10.0.txt')
ax.plot(XArr-Dx, AmpArr/amp,'k--', lw=5, markersize=10, label='R = 150 km' )
flength = np.append(flength, (XArr-Dx)[((AmpArr/amp)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp) ] - 3000.)
Amax = np.append(Amax, ((AmpArr/amp)[(XArr-Dx)> 3000.]).max())
flength2 = np.append(flength2, (XArr-Dx)[((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp*amp2/AmpArr2) ] - 3000.)
Aratio  = np.append(Aratio, ((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() )


XArr, AmpArr, amp, AmpArr2, amp2=get_data(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/homo_lens_10000km_2000km_R_200km_d_100km/field_data/Amp_10.0.txt')
ax.plot(XArr-Dx, AmpArr/amp,'b-', lw=5, markersize=10, label='R = 200 km' )
flength = np.append(flength, (XArr-Dx)[((AmpArr/amp)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp) ] - 3000.)
Amax = np.append(Amax, ((AmpArr/amp)[(XArr-Dx)> 3000.]).max())
flength2 = np.append(flength2, (XArr-Dx)[((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp*amp2/AmpArr2) ] - 3000.)
Aratio  = np.append(Aratio, ((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() )
# 
XArr, AmpArr, amp, AmpArr2, amp2=get_data(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/homo_lens_10000km_2000km_R_250km_d_100km/field_data/Amp_10.0.txt')
ax.plot(XArr-Dx, AmpArr/amp,'b--', lw=5, markersize=10, label='R = 250 km' )
flength = np.append(flength, (XArr-Dx)[((AmpArr/amp)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp) ] - 3000.)
Amax = np.append(Amax, ((AmpArr/amp)[(XArr-Dx)> 3000.]).max())
flength2 = np.append(flength2, (XArr-Dx)[((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp*amp2/AmpArr2) ] - 3000.)
Aratio  = np.append(Aratio, ((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() )

XArr, AmpArr, amp, AmpArr2, amp2=get_data(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/homo_lens_10000km_2000km_R_300km_d_100km/field_data/Amp_10.0.txt')
ax.plot(XArr-Dx, AmpArr/amp,'r-', lw=5, markersize=10, label='R = 300 km' )
flength = np.append(flength, (XArr-Dx)[((AmpArr/amp)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp) ] - 3000.)
Amax = np.append(Amax, ((AmpArr/amp)[(XArr-Dx)> 3000.]).max())
flength2 = np.append(flength2, (XArr-Dx)[((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp*amp2/AmpArr2) ] - 3000.)
Aratio  = np.append(Aratio, ((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() )

XArr, AmpArr, amp, AmpArr2, amp2=get_data(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/homo_lens_10000km_2000km_R_350km_d_100km/field_data/Amp_10.0.txt')
ax.plot(XArr-Dx, AmpArr/amp,'r--', lw=5, markersize=10, label='R = 350 km' )
flength = np.append(flength, (XArr-Dx)[((AmpArr/amp)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp) ] - 3000.)
Amax = np.append(Amax, ((AmpArr/amp)[(XArr-Dx)> 3000.]).max())
flength2 = np.append(flength2, (XArr-Dx)[((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp*amp2/AmpArr2) ] - 3000.)
Aratio  = np.append(Aratio, ((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() )

XArr, AmpArr, amp, AmpArr2, amp2=get_data(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/homo_lens_10000km_2000km_R_400km_d_100km/field_data/Amp_10.0.txt')
ax.plot(XArr-Dx, AmpArr/amp,'g-', lw=5, markersize=10, label='R = 400 km' )
flength = np.append(flength, (XArr-Dx)[((AmpArr/amp)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp) ] - 3000.)
Amax = np.append(Amax, ((AmpArr/amp)[(XArr-Dx)> 3000.]).max())
flength2 = np.append(flength2, (XArr-Dx)[((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() == (AmpArr/amp*amp2/AmpArr2) ] - 3000.)
Aratio  = np.append(Aratio, ((AmpArr/amp*amp2/AmpArr2)[(XArr-Dx)> 3000.]).max() )


# print  (XArr-Dx)[((AmpArr/amp)[(XArr-Dx)> 1000.]).max() == (AmpArr/amp) ] - 1000., ((AmpArr/amp)[(XArr-Dx)> 1000.]).max()
# XArr, AmpArr, amp, AmpArr2, amp2=get_data(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/cosine_anomaly_8000km_2000km_R_500km/field_data/Amp_10.0.txt')
# ax.plot(XArr-Dx, AmpArr/amp*amp2/AmpArr2,'m-', lw=3, markersize=10, label='R = 500 km' )

ax.plot(XArr-Dx, AmpArr2/amp2,'g--', lw=5, markersize=10, label='reference' );
plt.legend(loc='upper right', fontsize=20, numpoints = 1)
# y1=2900.
# y2=3100.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='k', alpha=0.5)
# y1=2800.
# y2=2900.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='b', alpha=0.5)
# y1=3100.
# y2=3200.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='b', alpha=0.5)
# y1=2700.
# y2=2800.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='r', alpha=0.5)
# y1=3200.
# y2=3300.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='r', alpha=0.5)
# y1=2600.
# y2=2700.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='g', alpha=0.5)
# y1=3300.
# y2=3400.
# ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='g', alpha=0.5)

ax.tick_params(axis='x', labelsize=25)
ax.tick_params(axis='y', labelsize=25)
# plt.yticks(np.arange(1., 3.2, 0.2))
plt.ylim([.2, 1.2])
# plt.xticks(np.arange(200., 7700., 500.))
plt.xticks(np.arange(500., 7500., 500.))
plt.xlim([500., 6000.])
# plt.ylim([(TArr-TArr2-1.5).min(), (TArr-TArr2).max()])
# plt.ylabel('Normalized amplitude', fontsize=30);
plt.ylabel('Amplitude ', fontsize=40);
# plt.xlabel('X position (km)', fontsize=30);
plt.xlabel('Distance (km)', fontsize=40);
# plt.title('Amplitude measurement with rectangle anomaly', fontsize=30);
# plt.title('Amplitude measurement with circular anomaly', fontsize=30);
# 
# 
# 
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
# RArr=np.array([200, 300,400])
RArr=np.array([100, 150, 200, 250, 300, 350, 400])
ax1.plot(RArr, Amax, 'r-o', lw=3, ms=20, label='Maximum amplitude' )
ax1.set_ylabel('Maximum amplitude', color='r', fontsize=40)
ax1.tick_params(axis='y', labelsize=35, color='r')
ax1.tick_params(axis='x', labelsize=35, color='k')
ax2.plot(RArr, Aratio, 'b-o', lw=3, ms=20, label='Maximum amplitude ratio' )
ax2.set_ylabel('Maximum amplitude ratio', fontsize=40, color='b');
ax1.set_xlabel('Radius of anomaly (km)', fontsize=40);
ax2.tick_params(axis='y', labelsize=35)
ax1.set_ylim(0.6, 1.0)
ax2.set_ylim(1.6, 3.0)

# # 
# # # plt.xticks(np.arange(200., 7700., 500.))
# # # plt.xlim([200., 7200.])
# # # plt.ylim([(TArr-TArr2-1.5).min(), (TArr-TArr2).max()])
# # # plt.ylabel('Normalized amplitude', fontsize=30);
# # plt.ylabel('Amplitude ratio', fontsize=40);
# # # plt.xlabel('X position (km)', fontsize=30);
# # plt.xlabel('Anomaly diameter (km)', fontsize=40);
# # 
fig, ax=plt.subplots()
# RArr=np.array([100, 200, 300,400])
plt.plot(RArr, flength2, 'k-o', lw=3, ms=20)
plt.ylabel('Focal length (km)', fontsize=40);
plt.xlabel('Radius of anomaly (km)', fontsize=40);
ax.tick_params(axis='x', labelsize=35)
ax.tick_params(axis='y', labelsize=35)
# 
# # plt.xticks(np.arange(200., 7700., 500.))
# # plt.xlim([200., 7200.])
# # plt.ylim([(TArr-TArr2-1.5).min(), (TArr-TArr2).max()])
# # plt.ylabel('Normalized amplitude', fontsize=30);
# plt.ylabel('Focal length (km)', fontsize=40);
# # plt.xlabel('X position (km)', fontsize=30);
# plt.xlabel('Anomaly diameter (km)', fontsize=40);
print Amax, flength
plt.show()

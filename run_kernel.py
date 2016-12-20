
import field2d_cartesian
import numpy as np
import matplotlib.animation as animation
import copy


kfield=field2d_cartesian.kernel_field(xmin=0, xmax=4000000, Nx=800, \
                    zmin=0, zmax=2000000, Nz=400, datadir='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_006_adj/OUTPUT_FILES', nproc=12)
# kfield.read_kernel_file()
# kfield.GetElementIndex()
# kfield.SaveElementIndex('/work3/leon/kernel_homo')
# kfield.LoadElementIndex('/lustre/janus_scratch/life9360/specfem2d_working_dir')
# kfield.get_kernel_value()
# kfield.writeASDF(outfname='/lustre/janus_scratch/life9360/specfem2d_working_dir/multipathing_2000km_4000km_adj_002/kernel_T.h5')
# kfield2=copy.deepcopy(kfield)
# kfield3=copy.deepcopy(kfield)
# kfield.readASDF(infname='kernel_circular.h5')
kfield.readASDF(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_006_adj/kernel_T.h5')
# kfield3.readASDF(infname='kernel_rectangle.h5')
# kfield3.readASDF(infname='kernel_rectangle.h5')
# kfield4=copy.deepcopy(kfield2)
# kfield.kernelvalue=kfield.kernelvalue-kfield2.kernelvalue
kfield.plot_kernel(vmin=-1e-8, vmax=1e-8)
# dd=kfield.get_dvalue( Xc=1000000, Zc=1000000, R=100000, vs=3000, dv = -0.1)
# kfield.writeASDF(outfname='kernel_rectangle.h5')

# 
# kfield.get_value(x=1500000)
# kfield2.get_value(x=1500000)
# 
# import matplotlib.pyplot as plt
# 
# fig, ax=plt.subplots()
# n=3; wavelength=30
# Fr=np.sqrt(n*wavelength*(1500.-200.)*(3000.-1500.)/(3000.-200.))
# Zarr=kfield.ZArr[:,0]
# plt.plot(Zarr/1000, kfield.value, 'k-', lw=3, markersize=10, label='circular')
# plt.plot(Zarr/1000, kfield2.value, 'r--',lw=3,  markersize=10, label='homo')
# ax.fill_betweenx(np.array([-5e-9, 2e-9]), 1000-Fr, 1000+Fr, facecolor='yellow', alpha=0.5)
# # plt.plot(Fr*np.ones(2)+1000, np.array([-5e-9, 2e-9]), 'r-',lw=3,  markersize=10, label='1st Fresnel')
# # plt.plot(1000-Fr*np.ones(2), np.array([-5e-9, 2e-9]), 'r-',lw=3,  markersize=10)
# 
# n=2; wavelength=30
# Fr=np.sqrt(n*wavelength*(1500.-200.)*(3000.-1500.)/(3000.-200.))
# ax.fill_betweenx(np.array([-5e-9, 2e-9]), 1000-Fr, 1000+Fr, facecolor='green', alpha=0.5)
# # plt.plot(Fr*np.ones(2)+1000, np.array([-5e-9, 2e-9]), 'g-',lw=3,  markersize=10, label='2nd Fresnel')
# # plt.plot(1000-Fr*np.ones(2), np.array([-5e-9, 2e-9]), 'g-',lw=3,  markersize=10)
# 
# n=1; wavelength=30
# Fr=np.sqrt(n*wavelength*(1500.-200.)*(3000.-1500.)/(3000.-200.))
# # plt.plot(Fr*np.ones(2)+1000, np.array([-5e-9, 2e-9]), 'b-',lw=3,  markersize=10, label='3rd Fresnel')
# # plt.plot(1000-Fr*np.ones(2), np.array([-5e-9, 2e-9]), 'b-',lw=3,  markersize=10)
# ax.fill_betweenx(np.array([-5e-9, 2e-9]), 1000-Fr, 1000+Fr, facecolor='blue', alpha=0.5)
# plt.legend(loc='lower right', fontsize=25)
# plt.show()
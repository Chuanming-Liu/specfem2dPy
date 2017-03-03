
import field2d_cartesian
import numpy as np
import matplotlib.animation as animation
import copy


kfield=field2d_cartesian.kernel_field(xmin=0, xmax=4000000, Nx=800, \
                    zmin=0, zmax=2000000, Nz=400, datadir='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_005_adj/OUTPUT_FILES', nproc=12)
# kfield.read_kernel_file()
# kfield.GetElementIndex()
# kfield.SaveElementIndex('/work3/leon/kernel_homo')
# kfield.LoadElementIndex('/lustre/janus_scratch/life9360/specfem2d_working_dir')
# kfield.get_kernel_value()
# kfield.writeASDF(outfname='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_005_adj/kernel_amp.h5')
# kfield2=copy.deepcopy(kfield)
# kfield3=copy.deepcopy(kfield)
# kfield.readASDF(infname='kernel_circular.h5')
kfield.readASDF(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_006_adj/kernel_amp.h5')
# kfield.readASDF(infname='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_006_adj/kernel_T.h5')
# kfield3.readASDF(infname='kernel_rectangle.h5')
# kfield3.readASDF(infname='kernel_rectangle.h5')
# kfield4=copy.deepcopy(kfield2)
# kfield.kernelvalue=kfield.kernelvalue-kfield2.kernelvalue
# kfield.plot_kernel(vmin=-1e-8, vmax=1e-8)
# kfield.plot_kernel(vmin=-2e-9, vmax=2e-9)
# dd=kfield.get_dvalue( Xc=1000000, Zc=1000000, R=100000, vs=3000, dv = -0.1)
# kfield.writeASDF(outfname='kernel_rectangle.h5')

# # 
kfield.get_value(x=1600000)
# kfield2.get_value(x=1500000)
# 
import matplotlib.pyplot as plt
# 
fig, ax=plt.subplots()
n=3; wavelength=30
Fr=np.sqrt(n*wavelength*(1500.-200.)*(3000.-1500.)/(3000.-200.))
Zarr=kfield.ZArr[:,0]
# plt.plot(Zarr/1000, kfield.value, 'k-', lw=3, markersize=10, label='circular')
# plt.plot(Zarr/1000, kfield.value*1e8, 'k-',lw=5,  markersize=10, label='1500')
# plt.ylabel(r"${\mathrm{10}^\mathrm{-8}}{\mathrm{s}}\cdot{\mathrm{m}^\mathrm{-2}}$", fontsize=35)

plt.plot(Zarr/1000, kfield.value*1e9, 'k-',lw=5,  markersize=10, label='1500')
plt.ylabel(r"${\mathrm{10}^\mathrm{-9}}{\mathrm{m}^\mathrm{-2}}$", fontsize=35)
unit='km'
plt.xlabel('z('+unit+')', fontsize=35)
plt.yticks(fontsize=30)
plt.xticks(fontsize=30)
        

plt.show()
#!/usr/bin/env python
import field2d_cartesian
import numpy as np
import matplotlib.animation as animation

# ## Wavefield Snapshots
Dx=1000.*100.
Dz=1000.*1000.
# WS=field2d_cartesian.WaveSnapshot(xmin=0, xmax=2000000, Nx=400, zmin=0, zmax=2000000, Nz=400, datadir=
#         '/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_D_570/OUTPUT_FILES', dn=50, nt=10000)
WS=field2d_cartesian.WaveSnapshot(xmin=0, xmax=4000000, Nx=800, zmin=0, zmax=2000000, Nz=400, datadir=
        '/lustre/janus_scratch/life9360/specfem2d_working_dir/multipathing_2000km_4000km/OUTPUT_FILES', dn=50, nt=25000)
# WS=field2d_cartesian.WaveSnapshot(xmin=0, xmax=1500000, Nx=300, zmin=0, zmax=1500000, Nz=300, datadir=
#          '/lustre/janus_scratch/life9360/LFMembrane_SH_D/OUTPUT_FILES', dn=50, nt=20000)
WS.ReadGridFile()

# WS.GetElementIndex()
# WS.SaveElementIndex('/lustre/janus_scratch/life9360/specfem2d_working_dir/multipathing_2000km_4000km')
WS.LoadElementIndex('/lustre/janus_scratch/life9360/specfem2d_working_dir/multipathing_2000km_4000km');
# WS.SaveElementIndex('/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_D/')
# WS.LoadElementIndex('/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_D')
WS.ReadSnapshots()
# im_ani=WS.PlotSnapshots()
# WS.writeASDF('field2d_d_001.h5')
# WS.ReadSingleSnap(1000)
# WS.PlotSingleSnap()
# WS.readASDF('field2d_d_001.h5')
# # WS.GetSingleSnap(1000)
# # WS.PlotSingleSnap()
# imani=WS.PlotSnapshots()
# # Writer = animation.writers['mencoder']
# imani.save('../specfem2d_4.mp4', fps=4, dpi=300)
# 
# # # 
# # # WS.LoadDSSnap()
# # # im_ani=WS.PlotDSSnaps()
# # # im_ani=WS.PlotDSSnaps(outfname='/projects/life9360/code/SEM/specfem2d/EXAMPLES/LFMembrane/homo.mp4')
# 

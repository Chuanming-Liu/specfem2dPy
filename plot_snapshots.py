#!/usr/bin/env python
import field2d_cartesian
import numpy as np
import matplotlib.animation as animation

# ## Wavefield Snapshots
Dx=1000.*100.
Dz=1000.*1000.
# WS=field2d_cartesian.WaveSnapshot(xmin=0, xmax=3000000+2*Dx, Nx=640, zmin=0, zmax=600000+2*Dz, Nz=520, datadir=
#         '/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_0.1_20/OUTPUT_FILES', dn=50, nt=10000)
WS=field2d_cartesian.WaveSnapshot(xmin=0, xmax=3000000+2*Dx, Nx=640, zmin=0, zmax=600000+2*Dz, Nz=520, datadir=
         '/lustre/janus_scratch/life9360/LFMembrane_SH_0.1_20/OUTPUT_FILES', dn=50, nt=10000)
# WS.ReadGridFile()
WS.GetElementIndex()
# WS.SaveElementIndex('/projects/life9360/code/SEM/specfem2d/EXAMPLES/LFMembrane_SH/OUTPUT_FILES');
# WS.LoadElementIndex('/projects/life9360/code/SEM/specfem2d/EXAMPLES/LFMembrane_SH/OUTPUT_FILES');
# WS.SaveElementIndex('/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_0.1_20/')
# WS.LoadElementIndex('/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_0.1_20')
# WS.ReadSnapshots();
# im_ani=WS.PlotSnapshots()
# WS.writeASDF('field2d_001.h5')
# WS.ReadSingleSnap(1000)
# WS.PlotSingleSnap()
# WS.readASDF('../field2d_001.h5')
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

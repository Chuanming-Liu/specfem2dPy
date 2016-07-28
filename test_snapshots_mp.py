
import field2d_cartesian
import numpy as np
import matplotlib.animation as animation

WS=field2d_cartesian.WaveSnapshot_mp(xmin=0, xmax=2000000, Nx=400, zmin=0, zmax=2000000, Nz=400, datadir=
        '/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_D_570/OUTPUT_FILES', dn=50, nt=16000, nproc=12)
# WS.ReadGridFile()
# WS.GetElementIndex()
# WS.SaveElementIndex('/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_D_570/OUTPUT_FILES')
# WS.LoadElementIndex('/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_D_570/OUTPUT_FILES')
# WS.ReadSingleSnap(2000)
# WS.ReadSnapshots()
# WS.writeASDF('../field2d_D_570.h5')
WS.readASDF('../field2d_D_570.h5')
WS.PlotSnapshots()

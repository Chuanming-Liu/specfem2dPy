
import field2d_cartesian
import numpy as np
import matplotlib.animation as animation

# WS=field2d_cartesian.WaveSnapshot_mp(xmin=0, xmax=2000000, Nx=400, zmin=0, zmax=2000000, Nz=400, datadir=
#         '/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_D_570/OUTPUT_FILES', dn=50, nt=16000, nproc=12)

WS=field2d_cartesian.WaveSnapshot_mp(xmin=0, xmax=4000000, Nx=800, zmin=0, zmax=2000000, Nz=400, datadir=
        '/lustre/janus_scratch/life9360/specfem2d_working_dir/2d_3d_staircase_basin_4000km_2000km_R_100km/OUTPUT_FILES', dn=50, nt=25000, nproc=12)
# WS.ReadGridFile()
# WS.GetElementIndex()
# WS.SaveElementIndex('/lustre/janus_scratch/life9360/specfem2d_working_dir/2d_3d_staircase_basin_4000km_2000km_R_100km/OUTPUT_FILES')
WS.LoadElementIndex('/lustre/janus_scratch/life9360/specfem2d_working_dir/2d_3d_staircase_basin_4000km_2000km_R_100km/OUTPUT_FILES')
# WS.ReadSingleSnap(6500)
WS.ReadSingleSnap(10300)
WS.PlotSingleSnap(factor=0.8)
# WS.ReadSnapshots()
# WS.writeASDF('/lustre/janus_scratch/life9360/specfem2d_working_dir/2d_3d_staircase_basin_4000km_2000km_R_100km/field2d_snap.h5')
# WS.readASDF('/lustre/janus_scratch/life9360/specfem2d_working_dir/2d_3d_staircase_basin_4000km_2000km/field2d_snap.h5')
# WS.PlotSnapshots()

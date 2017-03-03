#!/usr/bin/env python
import symdata
import matplotlib.pyplot as plt

dbase=symdata.specfem2dASDF('/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_004_adj/seismogram_adj.h5')
# 
dbase.readtxt('/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_004/DATA/STATIONS',
        datadir='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_004/OUTPUT_FILES')
dbase.AddEvent(x=200., z=1000.)


tr=dbase.get_trace(staid='MEM2D.300S0')
specTR=symdata.specfem2dtrace(tr.data, tr.stats)
specTR.get_adjoint_stf(outdir='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_004_adj', kerneltype='amp')
# specTR.get_adjoint_stf(outdir='/lustre/janus_scratch/life9360/specfem2d_working_dir/multipathing_2000km_4000km_adj_000', kerneltype='t', tmin=360, tmax=460)
#!/usr/bin/env python
import specfem2dPy
import numpy as np
import symData2d
import field2d_cartesian as field2d
import numpy.ma as ma
import matplotlib.pyplot as plt
SLst=specfem2dPy.StaLst();
SLst.ReadStaList('/projects/life9360/code/SEM/specfem2d/EXAMPLES/LFMembrane/DATA/STATIONS');
SDB=symData2d.Specfem2dDataBase(240000., 240000., SLst);

# 
# datadir='/projects/life9360/code/SEM/specfem2d/EXAMPLES/LFMembrane/OUTPUT_FILES'
# SDB.ReadtxtSeismograms(datadir=datadir);
# outdir='/lustre/janus_scratch/life9360/specfem2d_data/SAC_homo'
# SDB.SaveSeismograms(outdir);
# 
# # SDB.ReadSeismograms(datadir='/lustre/janus_scratch/life9360/specfem2d_data/SAC_homo')
# 
# outdir='/lustre/janus_scratch/life9360/specfem2d_data/aftan_homo';
# predV=np.array([ [0.1, 3.5], [0.5, 3.5], [1.0, 3.5], [5.0, 3.5], [20.0, 3.5]])
# 
# inftan=symData2d.InputFtanParam();
# inftan.setInParam(tmin=5.0, tmax=30.0, vmin=2.0, predV=predV);
# outdir='/lustre/janus_scratch/life9360/specfem2d_data/aftan_homo';
# SDB.aftanParallel(outdir=outdir, inftan=inftan);
# SDB.GetField2dFile(datadir='/lustre/janus_scratch/life9360/specfem2d_data/aftan_homo',
#                         outdir='/lustre/janus_scratch/life9360/specfem2d_data/field_homo', perLst=[10. ], outfmt='txt')
# 
dx=5
dy=5
ds=10
myfield=field2d.Field2d(Nx=480/dx, Ny=480/dx, dx=dx, dy=dy);
myfield.LoadFile('/lustre/janus_scratch/life9360/specfem2d_data/field_homo/TravelT.ph.10.0.txt')
myfield.natgridInterp();
# myfield.PlotField()
# myfield.Gradient();
# myfield.GetApparentV();
# myfield.PlotAppV()
myfield.GetDoT();
myfield.PlotDoT();
# myfield.LaplacianEqualXY()
# myfield.GetLplcCorrection(per=4.0)
# # myfield.Laplacian()
# # myfield.PlotLplcCo()

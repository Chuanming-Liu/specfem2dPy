#!/usr/bin/env python
import symData2d
import numpy as np
import field2d_cartesian as field2d
import matplotlib.pyplot as plt
import numpy.ma as ma
predV=np.array([ [0.1, 3.5], [0.5, 3.5], [1.0, 3.5], [5.0, 3.5], [10.0, 3.5]])
# predV=np.array([ [1.0, 2.1], [5.0, 2.5], [10.0, 2.9]])# 
dx=10
dy=10
# dx=5
# dy=5
# Database=symData2d.FatiandoDataBase(enx=500, enz=500, xmax=1000, zmax=1000, dx=dx, dz=dy, ds=1000);
Database=symData2d.FatiandoDataBase(enx=100, enz=500, xmax=601, zmax=601, dx=dx, dz=dy, ds=1000);
Database.ReadSeismograms(datadir='/lustre/janus_scratch/life9360/MEM2D_homo/SAC_003')
inftan=symData2d.InputFtanParam();
inftan.setInParam(tmin=2.0, tmax=10.0, vmin=2.0, predV=predV);
Database.aftanParallel(outdir='/lustre/janus_scratch/life9360/MEM2D_homo/DISP_aftan_003', inftan=inftan);
Database.GetField2dFile(datadir='/lustre/janus_scratch/life9360/MEM2D_homo/DISP_aftan_003',
                        outdir='/lustre/janus_scratch/life9360/MEM2D_homo/Field_aftan_003', perLst=[2., 3., 4., 6., 5., 7., 8., 9., 10. ], outfmt='txt')

# ds=10.
# myfield=field2d.Field2d(Nx=500/dx, Ny=500/dy, dx=ds, dy=ds);
# # myfield.LoadFile('/lustre/janus_scratch/life9360/MEM2D/Field_aftan/Amplitude.5.0.txt')
# # myfield.LoadFile('/lustre/janus_scratch/life9360/MEM2D/Field_freq/Amplitude.5.0.txt')
# myfield.LoadFile('/lustre/janus_scratch/life9360/MEM2D_slow/Field_aftan/Amplitude.4.0.txt')
# myfield.natgridInterp();
# # myfield.PlotField()
# 
# myfield.LaplacianEqualXY()
# myfield.GetLplcCorrection(per=4.0)
# # myfield.Laplacian()
# # myfield.PlotLplcCo()
# # 
# # myfield2=field2d.Field2d(Nx=500/dx, Ny=500/dy, dx=ds, dy=ds);
# # # myfield.LoadFile('/lustre/janus_scratch/life9360/MEM2D/Field_aftan/Amplitude.5.0.txt')
# # # myfield.LoadFile('/lustre/janus_scratch/life9360/MEM2D/Field_freq/Amplitude.5.0.txt')
# # myfield2.LoadFile('/lustre/janus_scratch/life9360/MEM2D_slow/Field_aftan/Amplitude.4.0.txt')
# # myfield2.natgridInterp();
# # # myfield2.PlotField()
# # 
# # myfield2.LaplacianEqualXY()
# # myfield2.GetLplcCorrection(per=4.0)
# # # myfield.Laplacian()
# # myfield2.PlotLplcCo()
# 
# 
ds=10.
# # dx=2
# # dy=2
Tfield=field2d.Field2d(Nx=601/dx, Ny=601/dy, dx=ds, dy=ds);
Tfield.LoadFile('/lustre/janus_scratch/life9360/MEM2D_homo/Field_aftan_003/TravelT.ph.5.0.txt')
Tfield.natgridInterp();
Tfield.Gradient();
Tfield.GetApparentV();
Tfield.PlotField()
# # # Tfield.GetCorV(myfield);
Tfield.PlotAppV()
plt.show()
# Tfield.LaplacianEqualXY();
# Tfield.PlotLaplacian()
# Tfield.PlotField()
# Tfield.PlotCorV()


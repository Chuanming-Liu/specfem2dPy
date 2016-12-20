#!/usr/bin/env python
import symdata



dbase=symdata.specfem2dASDF('/lustre/janus_scratch/life9360/specfem2d_working_dir/kernel_travelT_0.001/seismogram.h5')

dbase.readtxt('/lustre/janus_scratch/life9360/specfem2d_working_dir/kernel_travelT_0.001/DATA/STATIONS',
        datadir='/lustre/janus_scratch/life9360/specfem2d_working_dir/kernel_travelT_0.001/OUTPUT_FILES', factor=10)
dbase.AddEvent(x=200., z=1000.)
inftan=symdata.InputFtanParam()
inftan.ffact=5.
inftan.pmf=True
# inftan.vmin=2.5
# inftan.vmax=4.0
try:
    del dbase.auxiliary_data.DISPbasic1
    del dbase.auxiliary_data.DISPbasic2
    del dbase.auxiliary_data.DISPpmf1
    del dbase.auxiliary_data.DISPpmf2
except:
    pass
# 
# dbase.aftanMP(outdir='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_D_780/DATA/DISP', inftan=inftan, basic2=True,
#             pmf1=True, pmf2=True, tb=-11.9563813705)
dbase.aftan(inftan=inftan, basic2=True,
            pmf1=True, pmf2=True, tb=-11.9563813705)
# dbase.aftan(tb=-11.9563813705, outdir='/lustre/janus_scratch/life9360/LFMembrane_SH_0.1_20/DISP', inftan=inftan, basic2=True,
#             pmf1=True, pmf2=True)

try:
    del dbase.auxiliary_data.DISPpmf2interp
except:
    pass
dbase.InterpDisp(data_type='DISPbasic2')
# ndbase=dbase.SelectData(outfname='sw4synthetics001.h5', stafile='station_001.lst')

# del dbase.auxiliary_data.FieldDISPpmf2interp
dbase.GetField(outdir='/lustre/janus_scratch/life9360/specfem2d_working_dir/kernel_travelT_0.001/field_data', fieldtype='amp', data_type='DISPbasic2')
dbase.GetField(outdir='/lustre/janus_scratch/life9360/specfem2d_working_dir/kernel_travelT_0.001/field_data', fieldtype='Vgr',  data_type='DISPbasic2')
dbase.GetField(outdir='/lustre/janus_scratch/life9360/specfem2d_working_dir/kernel_travelT_0.001/field_data', fieldtype='Vph',  data_type='DISPbasic2')





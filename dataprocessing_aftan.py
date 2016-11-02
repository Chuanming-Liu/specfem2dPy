#!/usr/bin/env python
import symdata


dbase=symdata.specfem2dASDF('/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_005/seismogram.h5')

dbase.readtxt('/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_005/DATA/STATIONS',
        datadir='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_005/OUTPUT_FILES', factor=10)
dbase.AddEvent(x=200., z=1000.)
inftan=symdata.InputFtanParam()
inftan.ffact=5.
inftan.pmf=True
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
dbase.InterpDisp(data_type='DISPpmf2')
# ndbase=dbase.SelectData(outfname='sw4synthetics001.h5', stafile='station_001.lst')

# inv=dbase.Readsac('ak135_station.lst',datadir='/home/lili/sw4_working_dir/ak135_001', comptype='u')
# dbase.AddEvent(x=1000, y=1000, z=0)

# del dbase.auxiliary_data.FieldDISPpmf2interp
dbase.GetField(outdir='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_005/field_data', fieldtype='amp', data_type='DISPpmf2')
dbase.GetField(outdir='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_005/field_data', fieldtype='Vgr',  data_type='DISPpmf2')
dbase.GetField(outdir='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_healing_005/field_data', fieldtype='Vph',  data_type='DISPpmf2')





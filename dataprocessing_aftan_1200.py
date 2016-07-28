#!/usr/bin/env python
import symdata

# dbase=symdata.specfem2dASDF('synthetics_50.h5')
# dbase=symdata.specfem2dASDF('/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_D/specfem2d_570.h5')
dbase=symdata.specfem2dASDF('/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_D_1200/specfem2d_1200.h5')
# dbase=symdata.specfem2dASDF('/lustre/janus_scratch/life9360/LFMembrane_SH_0.1_20/specfem2d_mp.h5')
# dbase.readtxt('/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_D/OUTPUT_FILES/STATIONS',
#         datadir='/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_D/OUTPUT_FILES')

# dbase.readtxt('/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_D_1200/DATA/STATIONS',
#         datadir='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_D_1200/OUTPUT_FILES')
dbase.AddEvent(x=500., z=500.)
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
# dbase.aftanMP(outdir='/lustre/janus_scratch/life9360/specfem2d_working_dir/LFMembrane_SH_D_1200/DISP', inftan=inftan, basic2=True,
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
dbase.GetField(outdir='./D_1200km', fieldtype='amp', data_type='DISPpmf2')
dbase.GetField(outdir='./D_1200km', fieldtype='Vgr',  data_type='DISPpmf2')
dbase.GetField(outdir='./D_1200km', fieldtype='Vph',  data_type='DISPpmf2')





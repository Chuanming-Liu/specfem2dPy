#!/usr/bin/env python
import symdata

# dbase=symdata.specfem2dASDF('synthetics_50.h5')
dbase=symdata.specfem2dASDF('/lustre/janus_scratch/life9360/LFMembrane_SH_0.1_20/specfem2d.h5')
dbase.readtxt('STATIONS_z_300',
        datadir='/lustre/janus_scratch/life9360/LFMembrane_SH_0.1_20/OUTPUT_FILES')
dbase.AddEvent(x=1600., z=1300.)
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

# del dbase.auxiliary_data.FieldDISPbasic1interp
dbase.GetField(outdir='./ak135_4mem2d_20km_0.1', fieldtype='amp', data_type='DISPpmf2')
dbase.GetField(outdir='./ak135_4mem2d_20km_0.1', fieldtype='Vgr',  data_type='DISPpmf2')






#!/usr/bin/env python
import symdata

# dbase=symdata.specfem2dASDF('synthetics_50.h5')
dbase=symdata.specfem2dASDF('../specfem2d_0.1.h5')
# dbase.readtxt('station_downsample.lst',
#         datadir='/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH/OUTPUT_FILES')
# dbase.AddEvent(x=500., z=500.)
inftan=symdata.InputFtanParam()
inftan.ffact=1.
inftan.pmf=False
# try:
#     del dbase.auxiliary_data.DISPbasic1
#     del dbase.auxiliary_data.DISPbasic2
#     del dbase.auxiliary_data.DISPpmf1
#     del dbase.auxiliary_data.DISPpmf2
# except:
#     pass
# dbase.aftan( tb=-11.9563813705, inftan=inftan, basic2=False,
#             pmf1=False, pmf2=False)

# try:
#     del dbase.auxiliary_data.DISPbasic1interp
#     del dbase.auxiliary_data.DISPpmf2interp
# except:
#     pass
# dbase.InterpDisp(data_type='DISPbasic1')
# dbase.InterpDisp(data_type='DISPpmf2')
# dbase.GetField(outdir='./10km', fieldtype='amp', data_type='DISPbasic1')
# dbase.GetField(outdir='./10km', fieldtype='Vgr',  data_type='DISPbasic1', verbose=True)
# dbase.GetField(outdir='./10km', fieldtype='Vph',  data_type='DISPbasic1')
# dbase.GetField(outdir='./10km', fieldtype='amp', data_type='DISPbasic1')






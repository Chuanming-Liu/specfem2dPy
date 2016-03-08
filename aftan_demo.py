import obspy
import symData2d
# import ses3d_noisepy as snpy
import matplotlib.pyplot as plt
import numpy as np
f1='/lustre/janus_scratch/life9360/specfem2d_data/SAC_homo/MEM2D.35S5.SAC'
# f1='/projects/life9360/instaseis_seismogram/S001_INSTASEIS.LXZ.SAC'
# f1='/lustre/janus_scratch/life9360/EA_postprocessing_AGU/NKNT/EA.9100S4200..BXZ.SAC'
# f1='/Users/leon/MEM2d/SAC_homo/MEM2D_325_275.SAC'
st=obspy.read(f1)
tr=st[0]
predV=np.array([ [0.1, 3.5], [0.5, 3.5], [10.0, 3.5], [50.0, 3.5]])
tr1=symData2d.symtrace(tr.data, tr.stats)
# tr1.getSpec();
# # tr1.plotfreq()
# t1=tr1.stats.starttime
# tr1.trim(t1+10, t1+30)
# tr1.plot()
# # t = tr1.stats.starttime
# # tr1.trim(t , t + 1300.)  
# # tr2=tr1.GaussianFilter(1);
# 
tr1.aftan(piover4=-1, pmf=False, vmin=1.5, vmax=6.5, tmin=5., tmax=40.0, predV=predV)
# # tr1.aftan()
# tr1.getSNRTime()
tr1.plotftan(plotflag=0)
plt.show()
# tr1.plotTfr();

# perLst=np.array([0.5, 1, 3, 5 ])
# perLst=perLst[::-1]
# tr1.getAmpFreq(perLst=perLst);
# tr1.getSpec()
# tr1.plotfreq()
# for i in np.arange(40):
#     print (i-20.)/10.
#     tr1.aftan(piover4=(i-20.)/10., pmf=False, vmin=2.5, vmax=4.5, tmin=1.0, tmax=5.0, predV=predV)
#     tr1.plotftan(plotflag=0)
#     plt.show()



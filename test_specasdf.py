
import symdata

# dbase=symdata.sw4ASDF('/home/lili/sw4synthetics_ak135_EX_z_1km.h5')
dbase=symdata.specfem2dASDF('../specfem2d.h5')
# dbase.Readsac('station.lst', '../syndata_dir_000/sac_dir/', verbose=True)
dbase.PlotStreamsDistance()

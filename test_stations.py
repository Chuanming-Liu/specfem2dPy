import stations
Dx=1000.*100.
Dz=1000.*1000.
SLst=stations.StaLst();
SLst.ReadStaList('STATIONS')
nSLst=SLst.SelectStations(z=300000+Dz, x0=1600000, z0=1300000., dist=100000.)
# n2=nSLst.GetStation('/lustre/janus_scratch/life9360/LFMembrane_SH_0.1_20/OUTPUT_FILES')
nSLst.WriteStaList('STATIONS_z_300')
# n2.WriteStaList('STATIONS_z_300_t')
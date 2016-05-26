import numpy as np    
a = np.loadtxt('/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_0.1_20/wavefield_grid_for_dumps_000.txt')

ind = np.lexsort((a[:,0],a[:,1]))    

a[ind]

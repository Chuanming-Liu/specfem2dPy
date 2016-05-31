import numpy as np    
a = np.loadtxt('/home/lili/code/specfem2d/EXAMPLES/LFMembrane_SH_0.1_20/wavefield_grid_for_dumps_000.txt')
ind = np.lexsort((a[:,0],a[:,1]))    
b=a[ind]
print 'start'
# Ntotal=b.size
# for i in np.arange(Ntotal/2):
#     if i%10000==0:
#         print 'Step:', i, 'of', Ntotal
#     aaa=b[i]
#     cc=np.floor(aaa/2500)
#     dd=np.remainder(aaa, 2500)
    # print aaa
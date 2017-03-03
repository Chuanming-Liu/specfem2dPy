import numpy as np

n = 1./0.9
R=400.
# r1=300.; r2 = -300.
r1=R; r2 = -R
d1=50.; d2=50.
# d1=400; d2=400.
d = d1+d2
do=3000.
f=1./ ( (n-1)*( 1/r1 - 1/r2 + (n-1)*d/r1/r2/n) )

di = 1./(1/f - 1/do)

print f, di, do+di
import numpy as np
knots = np.array([-1.0, -0.6546536707079772, 0.0, 0.6546536707079772, 1.0])
#knots=(0.5+0.5*knots)
element=np.array([0, 120000, 240000, 480000])
A=np.ones([5,4])
AAA=(A*element).T
BBB=(A*element).T
CCC=(0.5+0.5*A.T*knots)*120000
DDD=AAA+CCC
EEE=np.tile(DDD, 4)
FFF=np.repeat(EEE,5)

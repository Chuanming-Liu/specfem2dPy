from sympy.solvers import solve
from sympy import Symbol
import numpy as np
f = Symbol('f')

flength = np.array([]); diArr=np.array([]); RArr=np.array([100, 150, 200, 250, 300, 350, 400])

d0=3000.
di=220.
print solve(1/d0*f**2 + f - di, f)
print 1./(1./d0+1./di)
flength=np.append(flength, 1./(1./d0+1./di))
diArr=np.append(diArr, di)
di=380.
print solve(1/d0*f**2 + f - di, f)
print 1./(1./d0+1./di)
flength=np.append(flength, 1./(1./d0+1./di))
diArr=np.append(diArr, di)

di=580.
print solve(1/d0*f**2 + f - di, f)
print 1./(1./d0+1./di)
flength=np.append(flength, 1./(1./d0+1./di))
diArr=np.append(diArr, di)

di=740.
print solve(1/d0*f**2 + f - di, f)
print 1./(1./d0+1./di)
flength=np.append(flength, 1./(1./d0+1./di))
diArr=np.append(diArr, di)

di=920.
print solve(1/d0*f**2 + f - di, f)
print 1./(1./d0+1./di)
flength=np.append(flength, 1./(1./d0+1./di))
diArr=np.append(diArr, di)

di=1100.
print solve(1/d0*f**2 + f - di, f)
print 1./(1./d0+1./di)
flength=np.append(flength, 1./(1./d0+1./di))
diArr=np.append(diArr, di)

di=1320.
print solve(1/d0*f**2 + f - di, f)
print 1./(1./d0+1./di)
flength=np.append(flength, 1./(1./d0+1./di))
diArr=np.append(diArr, di)

import matplotlib.pyplot as plt
from scipy import stats

# slope, intercept, r_value, p_value, std_err = stats.linregress(RArr,flength)
# slope, intercept, r_value, p_value, std_err = stats.linregress(RArr,flength)

fig, ax = plt.subplots()
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(RArr,flength)
ax.plot(RArr, slope1 *RArr + intercept1, 'k-', lw=5, label='')

plt.ylabel('Focal length (km)', fontsize=40);
plt.xlabel('Radius of anomaly (km)', fontsize=40);
ax.plot(RArr, flength, 'o',  ms=20, label='Observed focal length')
ax.tick_params(axis='y', labelsize=35, color='k')
ax.tick_params(axis='x', labelsize=35, color='k')
print slope1, intercept1, r_value1, p_value1, std_err1

fig, ax = plt.subplots()
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(RArr,diArr)
ax.plot(RArr, slope2 *RArr + intercept2, 'k-', lw=5)
ax.plot(RArr, diArr, 'bo', ms=20)
ax.tick_params(axis='y', labelsize=35, color='k')
ax.tick_params(axis='x', labelsize=35, color='k')
plt.ylabel('Image distance (km)', fontsize=40);
plt.xlabel('Radius of anomaly (km)', fontsize=40);
print slope2, intercept2, r_value2, p_value2, std_err2
plt.show()


# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(RArr, diArr, 'ro', lw=3, ms=20, label='Image distance' )
# ax1.set_ylabel('Image distance (km)', color='r', fontsize=40)
# ax1.tick_params(axis='y', labelsize=35, color='r')
# ax1.tick_params(axis='x', labelsize=35, color='k')
# ax2.plot(RArr, flength, 'bo', lw=3, ms=20, label='Focal length (km)' )
# ax.plot(RArr, slope1 *RArr + intercept1, 'k-', lw=5)
# ax2.set_ylabel('Focal length (km)', fontsize=40, color='b');
# ax1.set_xlabel('Radius of anomaly (km)', fontsize=40);
# ax2.tick_params(axis='y', labelsize=35)
# # ax1.set_ylim(0.6, 1.0)
# # ax2.set_ylim(1.6, 3.0)
# plt.show()

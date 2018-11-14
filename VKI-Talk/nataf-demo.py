from equadratures import *
import numpy as np
import matplotlib.pyplot as plt


# Quick code that shows if you have a linear model, a Nataf transform
# may require a quadratic or higher to approximate output values

def fun(x):
    return 0.3*x[0] - 0.7*x[1]

myParam1 = Parameter(distribution='truncated-gaussian', shape_parameter_A=3.5, shape_parameter_B=2, lower=0., upper=6.0, order=1)
myParam2 = Parameter(distribution='truncated-gaussian', shape_parameter_A=-5., shape_parameter_B=0.1, lower=-6., upper=-2., order=1)

z, pdf = myParam1.getPDF()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.plot(z, pdf, '-', c='crimson', lw=4)
adjust_spines(ax, ['left', 'bottom'])
ax.set_axisbelow(True)
plt.xlabel('$\zeta_1$', fontsize=15)
plt.ylabel('PDF', fontsize=15)
plt.savefig('Figure_Nataf_Dist_1.png', dpi=200, bbox_inches='tight')
plt.show()

z, pdf = myParam2.getPDF()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.plot(z, pdf, '-', c='navy', lw=4)
adjust_spines(ax, ['left', 'bottom'])
ax.set_axisbelow(True)
plt.xlabel('$\zeta_2$', fontsize=15)
plt.ylabel('PDF', fontsize=15)
plt.savefig('Figure_Nataf_Dist_2.png', dpi=200, bbox_inches='tight')
plt.show()

# Define correlation matrix!
R = np.eye(2)
R[0, 1] = 0.7
R[1, 0] = 0.7


u1 = Parameter(distribution='normal', shape_parameter_A=0.0, shape_parameter_B=1.0, order=3)
myNataf = Nataf([myParam1, myParam2], R)

# For Monte-Carlo!
#samples_mc = myNataf.getCorrelatedSamples(N=10000)
#f_mc = evalfunction(samples_mc, fun)

# For Polynomials!
myBasis = Basis('Tensor grid')
myPoly = Polyint([u1, u1], myBasis)
samples_p =  myPoly.quadraturePoints
samples_corr_p = myNataf.U2C(samples_p)
f_p = evalfunction(samples_corr_p, fun)

myPoly.computeCoefficients(f_p)
myStats = myPoly.getStatistics()

print '----POLYNOMIALS-----'
print myStats.mean, myStats.variance, myStats.skewness
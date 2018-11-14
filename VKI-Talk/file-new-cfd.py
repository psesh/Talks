import numpy as np
from equadratures import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib
params = {
          'axes.facecolor':'white',
         }
matplotlib.rcParams.update(params)


zeta_1 = Parameter(distribution='uniform', order=4, lower= -2.0, upper=2.0)
zeta_2 = Parameter(distribution='uniform', order=4, lower=-2.0, upper=2.5)

def fun(x):
        return x[0] - 2*x[1]**2*x[0] - 3*x[0]**4 + 88.3


#myBasis = Basis('total order', [zeta_1.order, zeta_2.order])
myBasis = Basis('hyperbolic basis', [zeta_1.order, zeta_2.order], q=0.7)
myPoly = Polylsq([zeta_1, zeta_2], myBasis, mesh='tensor', optimization='greedy-qr', oversampling=1.0)
myPoly.computeCoefficients(fun)
coefficients = myPoly.coefficients
Az = myPoly.Az
m, n = Az.shape
P = np.linalg.pinv(Az)
Sigma = np.eye((m)) * 0.3
Sigma[0,0] = 0.4
uncertainties = np.sqrt(np.diag(Sigma))*1.96
cov_in_coefficients = np.dot( np.dot(P, Sigma), P.T)
std_in_coefficients =  np.reshape( np.sqrt( np.diag(cov_in_coefficients) ) * 1.96, (n, 1) )
quadraturePoints = myPoly.quadraturePoints
myPolyLow = Poly([zeta_1, zeta_2], myBasis)
myPolyLow.__setCoefficients__(coefficients - std_in_coefficients)
myPolyHigh = Poly([zeta_1, zeta_2], myBasis)
myPolyHigh.__setCoefficients__(coefficients + std_in_coefficients)

# For plotting!
N = 20
z1 = np.linspace(zeta_1.lower, zeta_1.upper, N)
z2 = np.linspace(zeta_2.lower, zeta_2.upper, N)
[Z1, Z2] = np.meshgrid(z1, z2)
Z1_vec = np.reshape(Z1, (N*N, 1))
Z2_vec = np.reshape(Z2, (N*N, 1))
samples = np.hstack([Z1_vec, Z2_vec])

PolyApprox = np.reshape( myPoly.evaluatePolyFit( samples ) , (N, N) ) 
PolyApproxLow = np.reshape( myPolyLow.evaluatePolyFit( samples ) , (N, N) ) 
PolyApproxHigh = np.reshape( myPolyHigh.evaluatePolyFit( samples ) , (N, N) ) 

#PolyApprox = myPoly.evaluatePolyFit( samples )
#PolyApprox = np.reshape(PolyApprox, (N, N))

# Response surfaces
zval = myPoly.functionEvaluations
print len(zval)
print len(uncertainties)
print len(quadraturePoints[:,0])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(Z1, Z2, PolyApprox, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, alpha=0.8)
wire = ax.plot_surface(Z1, Z2, PolyApproxLow, color='gray', rstride=1, cstride=1, alpha=0.3, linewidth=0.2)
wire = ax.plot_surface(Z1, Z2, PolyApproxHigh, color='gray', rstride=1, cstride=1, alpha=0.3, linewidth=0.2)
ax.scatter(quadraturePoints[:,0], quadraturePoints[:,1], myPoly.functionEvaluations, c='gold', s=80, edgecolor='black')
for i in range(0, m):
    xvec = [quadraturePoints[i,0] , quadraturePoints[i,0]] 
    yvec = [quadraturePoints[i,1] , quadraturePoints[i,1]] 
    zvec = [ float(zval[i]- uncertainties[i]), float(zval[i] + uncertainties[i])] 
    ax.plot( xvec, yvec, zvec, color='black' )

ax.set_xlabel('$\zeta_1$', fontsize=15)
ax.set_ylabel('$\zeta_2$', fontsize=15)
ax.set_zlabel('$f(\zeta_1, \zeta_2)$', fontsize=15)
ax.view_init(7, 72)
frame1 = plt.gca()
ax.set_yticks([])
ax.set_xticks([])
ax.set_zticks([])
ax.grid(True)
frame1.axes.yaxis.set_ticklabels([])
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.zaxis.set_ticklabels([])
plt.savefig('Fig_007.png', dpi=200, bbox_inches='tight')
plt.show()
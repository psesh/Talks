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


myBasis = Basis('hyperbolic basis', q=0.7)
myPoly = Polylsq([zeta_1, zeta_2], myBasis, mesh='tensor', optimization='greedy-qr', oversampling=1.0)
myPoly.computeCoefficients(fun)

myBasis2 = Basis('tensor grid')
myPoly2 = Polyint([zeta_1, zeta_2], myBasis2)
myPoly2.computeCoefficients(fun)

effectiveQuadraturePoints = myPoly.quadraturePoints
tensorQuadraturePoints = myPoly2.quadraturePoints

# For plotting!
N = 20
z1 = np.linspace(zeta_1.lower, zeta_1.upper, N)
z2 = np.linspace(zeta_2.lower, zeta_2.upper, N)
[Z1, Z2] = np.meshgrid(z1, z2)
Z1_vec = np.reshape(Z1, (N*N, 1))
Z2_vec = np.reshape(Z2, (N*N, 1))
samples = np.hstack([Z1_vec, Z2_vec])

PolyApprox = np.reshape( myPoly.evaluatePolyFit( samples ) , (N, N) ) 
PolyApprox2 = np.reshape( myPoly2.evaluatePolyFit( samples ) , (N, N) ) 


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(Z1, Z2, PolyApprox, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, alpha=0.8)
ax.scatter(effectiveQuadraturePoints[:,0], effectiveQuadraturePoints[:,1], myPoly.functionEvaluations, c='gold', s=80, edgecolor='black')
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
plt.savefig('Fig_007_A.png', dpi=200, bbox_inches='tight')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(Z1, Z2, PolyApprox2, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, alpha=0.8)
ax.scatter(tensorQuadraturePoints[:,0], tensorQuadraturePoints[:,1], evalfunction(tensorQuadraturePoints, fun), c='gold', s=80, edgecolor='black')
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
plt.savefig('Fig_007_B.png', dpi=200, bbox_inches='tight')
plt.show()

print myPoly.basis.cardinality
print myPoly2.basis.cardinality


# Original!
A = myPoly.A
G = np.dot(A.T, A)
mm, nn = A.shape
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.imshow(G, cmap=cm.jet)
plt.colorbar()
frame1 = plt.gca()
ax.set_xticks([])
frame1.axes.xaxis.set_ticklabels([])
ticks = []
for i in range(0, myBasis.cardinality):
    ticks.append(str(myBasis.elements[i,:]))
plt.yticks(np.arange(0, myBasis.cardinality, 2), ticks[0: myBasis.cardinality:2], fontsize=8)
#plt.yticks(np.arange(0, small_tensor.cardinality, 3), ticks[0:small_tensor.cardinality:3], fontsize=10)
ax.grid(False)
titlestring = r'$\mathbf{G} = \mathbf{A}^T \mathbf{A}$ where $\mathbf{A}\in \mathbb{R}^{%s \times %s}$'%(mm, nn)
plt.title(titlestring,fontsize=16)
plt.savefig('Gram_Orig_007A.png' , dpi=200, bbox_inches='tight')
plt.show()


# Original!
Az = myPoly.Az
Gz = np.dot(Az.T , Az)
mm, nn = Az.shape
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.imshow(Gz, cmap=cm.jet)
plt.colorbar()
frame1 = plt.gca()
ax.set_xticks([])
frame1.axes.xaxis.set_ticklabels([])
ticks = []
for i in range(0, myBasis.cardinality):
    ticks.append(str(myBasis.elements[i,:]))
plt.yticks(np.arange(0, myPoly.basis.cardinality, 2), ticks[0: myPoly.basis.cardinality:2], fontsize=8)
ax.grid(False)
titlestring = r'$\mathbf{Gz} = \mathbf{Az}^T \mathbf{Az}$ where $\mathbf{Az}\in \mathbb{R}^{%s \times %s}$'%(mm, nn)
plt.title(titlestring,fontsize=16)
plt.savefig('Gram_Mod_007B.png' , dpi=200, bbox_inches='tight')
plt.show()

print np.linalg.cond(G)
print np.linalg.cond(Gz)

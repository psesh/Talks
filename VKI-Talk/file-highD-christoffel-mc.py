import numpy as np
from equadratures import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib
from scipy.stats import arcsine
from scipy.linalg import qr, svd, lu
params = {
          'axes.facecolor':'white',
         }
matplotlib.rcParams.update(params)



zeta_1 = Parameter(distribution='uniform', order=3, lower= -1.0, upper=1.0)
zeta_2 = Parameter(distribution='uniform', order=3, lower=-1.0, upper=1.0)
zeta_3 = Parameter(distribution='uniform', order=3, lower=-1.0, upper=1.0)

myBasis = Basis('Total order', orders=[4, 4, 4])
myPoly = Poly([zeta_1, zeta_2, zeta_3], myBasis)

no_of_points = 1000

# Christoffel!
p = arcsine.rvs(size=(no_of_points, 3))*2.0 - 1
P = myPoly.getPolynomial(p)
w =  (no_of_points * myBasis.cardinality)/np.sum( P**2 , 0)
w = w * 1.0/np.sum(w)
W = np.diag(np.sqrt(w))
A = np.dot(W, P.T)
G = np.dot(A.T, A)
mm, nn = A.shape
print np.linalg.cond(G)

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
plt.yticks(np.arange(0, myBasis.cardinality, 4), ticks[0: myBasis.cardinality:4], fontsize=9)
#plt.yticks(np.arange(0, small_tensor.cardinality, 3), ticks[0:small_tensor.cardinality:3], fontsize=10)
ax.grid(False)
titlestring = r'$\mathbf{G} = \mathbf{A}^T \mathbf{A}$ where $\mathbf{A}\in \mathbb{R}^{%s \times %s}$'%(mm, nn)
plt.title(titlestring,fontsize=16)
plt.savefig('Christoffel_image.png' , dpi=200, bbox_inches='tight')
plt.show()



# Standard MC!
p = np.random.rand(no_of_points, 3 )*2.0 - 1 #arcsine.rvs(size=(no_of_points, 5))*2.0 - 1
P = myPoly.getPolynomial(p)
w = (no_of_points * myBasis.cardinality)/np.sum( P**2 , 0)


w = w * 1.0/np.sum(w)
W = np.diag(np.sqrt(w))
A = np.dot(W, P.T)
G = np.dot(A.T, A)
mm, nn = A.shape
print np.linalg.cond(G)

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
plt.yticks(np.arange(0, myBasis.cardinality, 4), ticks[0: myBasis.cardinality:4], fontsize=9)
#plt.yticks(np.arange(0, small_tensor.cardinality, 3), ticks[0:small_tensor.cardinality:3], fontsize=10)
ax.grid(False)
titlestring = r'$\mathbf{G} = \mathbf{A}^T \mathbf{A}$ where $\mathbf{A}\in \mathbb{R}^{%s \times %s}$'%(mm, nn)
plt.title(titlestring,fontsize=16)
plt.savefig('MC_image.png' , dpi=200, bbox_inches='tight')
plt.show()



#zhat, L, ztilde, Utilde = maxdet(A, nn)
#z = binary2indices(zhat)
"""
__, __, pvec = qr(A.T, pivoting=True)
z = pvec[0:nn]

points_subsampled = p[z]
wts_orig_normalized =  w[z] / np.sum(w[z]) # if we pick a subset of the weights, they should add up to 1.!
Pz = myPoly.getPolynomial(points_subsampled)
Wz = np.mat(np.diag( np.sqrt(wts_orig_normalized) ) )
Az =  Wz * Pz.T
Gz = np.dot(Az.T, Az)
print np.linalg.cond(Gz)
print Gz.shape 
"""



"""
# Original!
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
plt.yticks(np.arange(0, myPoly.basis.cardinality, 4), ticks[0: myPoly.basis.cardinality:4], fontsize=9)
ax.grid(False)
titlestring = r'$\mathbf{Gz} = \mathbf{Az}^T \mathbf{Az}$ where $\mathbf{Az}\in \mathbb{R}^{%s \times %s}$'%(mm, nn)
plt.title(titlestring,fontsize=16)
plt.savefig('Gram_Mod_007C.png' , dpi=200, bbox_inches='tight')
plt.show()

print np.linalg.cond(G)
print np.linalg.cond(Gz)
"""

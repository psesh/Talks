from equadratures import *
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm
from scipy.stats import arcsine
from scipy.linalg import qr, svd, lu
from equadratures.convex import maxdet, binary2indices

myParam = Parameter(distribution='uniform', lower=-1., upper=1., order=10)
myBasis = Basis('Univariate', orders=[10])
myBasisNew = Basis('Univariate', orders=[500])
myPoly = Polyint([myParam], myBasis)

p2, w2 = myParam._getLocalQuadrature()
no_of_points = 500
p = arcsine.rvs(size=no_of_points)*2.0 - 1
Polymat = myPoly.getPolynomial(p)
w =  1.0/np.sum( Polymat**2 , 0)
w = w * 1.0/np.sum(w)
P = myPoly.getPolynomial(p)
W = np.diag(np.sqrt(w))
A = np.dot(W, P.T)
G = np.dot(A.T, A)
mm, nn = A.shape
print np.linalg.cond(G)


#zhat, L, ztilde, Utilde = maxdet(A, nn)
#z = binary2indices(zhat)

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

# Original!
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
plt.savefig('Gram_Orig_A1.png' , dpi=200, bbox_inches='tight')
plt.show()


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
plt.yticks(np.arange(0, myBasis.cardinality, 2), ticks[0: myBasis.cardinality:2], fontsize=8)
ax.grid(False)
titlestring = r'$\mathbf{Gz} = \mathbf{Az}^T \mathbf{Az}$ where $\mathbf{Az}\in \mathbb{R}^{%s \times %s}$'%(mm, nn)
plt.title(titlestring,fontsize=16)
plt.savefig('Gram_Mod_A1.png' , dpi=200, bbox_inches='tight')
plt.show()


print p2.shape 

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.scatter(p, np.zeros((1,500))+0.1,  marker='o', s=10,  color='navy', linewidth=1., label='Arcsine samples')
plt.scatter(points_subsampled, np.zeros((11, 1)),  marker='o', s=80,  color='crimson',linewidth=2., label='QR-pivoted subsamples')
plt.scatter(p2, np.zeros((11, 1))-0.1,  marker='o', s=80,  color='gold',linewidth=2., label='Gauss-Legendre', edgecolor='black')
adjust_spines(ax, ['left', 'bottom'])
plt.xlabel('Points, $\omega_i$', fontsize=13)
#titlestring = 'Sampling, $m$=%s'%(myParam.order + 1)
#plt.title(titlestring,fontsize=16)
frame1 = plt.gca()
ax.set_yticks([])
frame1.axes.yaxis.set_ticklabels([])
plt.xlim([-1.05, 1.05])
plt.ylim([-0.2, 0.2])
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0., fontsize=12)
plt.savefig('Points_new_gram.png', dpi=200, bbox_inches='tight', pad_inches=0.1)
plt.show()

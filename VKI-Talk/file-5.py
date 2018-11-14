from equadratures import *
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm
from scipy.stats import arcsine

myParam = Parameter(distribution='uniform', lower=-1., upper=1., order=5)
myBasis = Basis('Univariate', orders=[50])
myBasisNew = Basis('Univariate', orders=[500])
myPoly = Polyint([myParam], myBasis)

p2, w2 = myParam._getLocalQuadrature()
no_of_points = 500
#p = arcsine.rvs(size=no_of_points)*2.0 - 1
p = np.random.rand(no_of_points, 1)*2 - 1.0

Polymat = myPoly.getPolynomial(p)
vv, zz = Polymat.shape

w =  1.0/np.sum( Polymat**2 , 0)
w = w * 1.0/np.sum(w)



P = myPoly.getPolynomial(p)
W = np.diag(np.sqrt(w))
A = np.dot(W, P.T)
print A.shape
G = np.dot(A.T, A)
print G.shape

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
plt.yticks(np.arange(0, myBasis.cardinality, 2), ticks[0: myBasis.cardinality:2], fontsize=8)
#plt.yticks(np.arange(0, small_tensor.cardinality, 3), ticks[0:small_tensor.cardinality:3], fontsize=10)
ax.grid(False)
titlestring = r'$\mathbf{G} = \mathbf{A}^T \mathbf{A}$ where $\mathbf{A}\in \mathbb{R}^{%s \times %s}$'%(mm, nn)
plt.title(titlestring,fontsize=16)
plt.savefig('Gram_B.png' , dpi=200, bbox_inches='tight')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
n, bins, patches = plt.hist(p, 20, density=True, edgecolor='white', linewidth=1.2, facecolor='crimson')
plt.grid(True)
adjust_spines(ax, ['left', 'bottom'])
titlestring = 'Sampling distribution'
plt.title(titlestring,fontsize=16)
plt.xlim([-1.05, 1.05])
plt.xlabel('Points, $\zeta_i$', fontsize=13)
plt.savefig('Distribution_B.png', dpi=200, bbox_inches='tight', pad_inches=0.1)
plt.show()

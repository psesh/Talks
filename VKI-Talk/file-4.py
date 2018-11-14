from equadratures import *
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm

myParam = Parameter(distribution='uniform', lower=-1., upper=1., order=10)
myBasis = Basis('Univariate', orders=[20])
myPoly = Polyint([myParam], myBasis)


p, w = myParam._getLocalQuadrature()

P = myPoly.getPolynomial(p)
W = np.diag(np.sqrt(w))
A = np.dot(W, P.T)
print A.shape
G = np.dot(A.T, A)
print G.shape


mm, nn = A.shape
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.imshow(np.log10(np.abs(G)), cmap=cm.jet, vmin=-13., vmax=0.0)
plt.colorbar()
frame1 = plt.gca()
ax.set_xticks([])
frame1.axes.xaxis.set_ticklabels([])
ticks = []
for i in range(0, myBasis.cardinality):
    ticks.append(str(myBasis.elements[i,:]))
plt.yticks(np.arange(myBasis.cardinality), ticks, fontsize=8)
ax.grid(False)
titlestring = r'$\mathbf{G} = \mathbf{A}^T \mathbf{A}$ where $\mathbf{A}\in \mathbb{R}^{%s \times %s}$'%(mm, nn)
plt.title(titlestring,fontsize=16)
plt.savefig('Gram_matrix_D.png' , dpi=200, bbox_inches='tight')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.scatter(p, w,  marker='o', s=100,  color='steelblue',linewidth=2.)
adjust_spines(ax, ['left', 'bottom'])
plt.xlabel('Weights, $\zeta_i$', fontsize=13)
plt.ylabel('Points, $\omega_i$', fontsize=13)
titlestring = 'Gauss-Legendre quadrature, $m$=%s'%(myParam.order + 1)
plt.title(titlestring,fontsize=16)
plt.xlim([-1.05, 1.05])
plt.ylim([0.0, 0.15])
plt.savefig('Points_D.png', dpi=200, bbox_inches='tight', pad_inches=0.1)
plt.show()

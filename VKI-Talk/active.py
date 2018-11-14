import numpy as np
import scipy
from equadratures import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib
params = {
          'axes.facecolor':'white',
         }
matplotlib.rcParams.update(params)

X = np.loadtxt('../../private-data/X.dat')
fX = np.loadtxt('../../private-data/Y.dat')

n=25
p=2


print X.shape

listing=[]
for i in range(0,n):
    listing.append(p)
basis=Basis('Total order',listing)

params=[]
P=Parameter(order=p,distribution='uniform',lower=-1,upper=1)
for i in range(0,n):
    params.append(P)

Polybasis= Polyreg(params, basis, training_inputs=X, training_outputs=fX)
[eigs, U]=computeActiveSubspaces(Polybasis, samples=None)


active1 = np.dot( X , U[:,0] )
active2 = np.dot( X , U[:,1] )

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
cax = ax.scatter(active1, active2, fX, c=fX, marker='o', s=50, edgecolor='black')
cbar = plt.colorbar(cax, extend='neither', spacing='proportional', orientation='vertical', shrink=0.8)
#ax.set_xlabel('Active variable 1, $\mathbf{XU_{1}}$')
#ax.set_ylabel('Non-dimensional efficiency, $\eta$')
plt.xlim([-2.0, 2.0])
ax.set_xlabel('$v_{1}^T \zeta_1$', fontsize=15)
ax.set_ylabel('$v_{2}^T \zeta_2$', fontsize=15)
ax.set_zlabel('$f(\zeta)$', fontsize=15)
ax.view_init(20, 104)
frame1 = plt.gca()
ax.set_yticks([])
ax.set_xticks([])
ax.set_zticks([])
ax.grid(True)
frame1.axes.yaxis.set_ticklabels([])
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.zaxis.set_ticklabels([])
plt.savefig('active_1.png' , dpi=200, bbox_inches='tight')
plt.show()

PolySobol = Polybasis.getStatistics(max_sobol_order=1)
data = PolySobol.getSobol(1)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for key, value in data.items():
    plt.bar(key[0], value, color='navy')
adjust_spines(ax, ['left', 'bottom'])
plt.xlabel("Parameters")
plt.ylabel("First order Sobol' indices")
plt.savefig('active_2.png' , dpi=200, bbox_inches='tight')
plt.show()

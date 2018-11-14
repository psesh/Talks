from equadratures import *
import numpy as np
import matplotlib.pyplot as plt


def fun(x):
    return np.exp(3*x[0] + x[1])


x = Parameter(distribution='uniform', lower=-1., upper=1., order=35)
tensor = Basis('Tensor grid')
Pol = Polyint([x,x], tensor)
Pol.computeCoefficients(fun)

x, y, z, max_order = twoDgrid(Pol.coefficients, Pol.multi_index)
G = np.log10(np.abs(z))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
cax = plt.scatter(x, y, s=30, marker='o', c=G, cmap='jet', alpha=1.0, vmin=-16.0, vmax=0.)
plt.xlim(-0.5, max_order)
plt.ylim(-0.5, max_order)
adjust_spines(ax, ['left', 'bottom'])
ax.set_axisbelow(True)
plt.xlabel('$i_1$', fontsize=13)
plt.ylabel('$i_2$', fontsize=13)
cbar = plt.colorbar(extend='neither', spacing='proportional',
                orientation='vertical', shrink=0.8, format="%.0f")
cbar.ax.tick_params(labelsize=13)
plt.savefig('Pseudo_1.png',   dpi=300, bbox_inches='tight')


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.scatter(Pol.quadraturePoints[:,0], Pol.quadraturePoints[:,1] , marker='o', s=2, color='tomato')
adjust_spines(ax, ['left', 'bottom'])
plt.xlabel('$\zeta_1$', fontsize=13)
plt.ylabel('$\zeta_2$', fontsize=13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.savefig('FigureB.png' , dpi=200, bbox_inches='tight', pad_inches=0.1)


del Pol
x = Parameter(distribution='uniform', lower=-1., upper=1., order=4)
sparse = Basis('Sparse grid', level=6, growth_rule='exponential')
Pol = Polyint([x,x], sparse)
Pol.computeCoefficients(fun)



x, y, z, max_order = twoDgrid(Pol.coefficients, Pol.multi_index)
G = np.log10(np.abs(z))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
cax = plt.scatter(x, y, s=30, marker='o', c=G, cmap='jet', alpha=1.0, vmin=-16.0, vmax=0.)
plt.xlim(-0.5, max_order)
plt.ylim(-0.5, max_order)
adjust_spines(ax, ['left', 'bottom'])
ax.set_axisbelow(True)
plt.xlabel('$i_1$', fontsize=13)
plt.ylabel('$i_2$', fontsize=13)
cbar = plt.colorbar(extend='neither', spacing='proportional',
                orientation='vertical', shrink=0.8, format="%.0f")
cbar.ax.tick_params(labelsize=13)
plt.savefig('Pseudo_2.png',   dpi=300, bbox_inches='tight')

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.scatter(Pol.quadraturePoints[:,0], Pol.quadraturePoints[:,1] , marker='o', s=2, color='tomato')
adjust_spines(ax, ['left', 'bottom'])
plt.xlabel('$\zeta_1$', fontsize=13)
plt.ylabel('$\zeta_2$', fontsize=13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.savefig('FigureD.png' , dpi=200, bbox_inches='tight', pad_inches=0.1)
print len(np.unique(Pol.quadraturePoints, axis=0))


del Pol
x = Parameter(distribution='uniform', lower=-1., upper=1., order=4)
sparse = Basis('Sparse grid', level=12, growth_rule='linear')
Pol = Polyint([x,x], sparse)
Pol.computeCoefficients(fun)
x, y, z, max_order = twoDgrid(Pol.coefficients, Pol.multi_index)
G = np.log10(np.abs(z))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
cax = plt.scatter(x, y, s=200, marker='o', c=G, cmap='jet', alpha=1.0, vmin=-16.0, vmax=0.)
plt.xlim(-0.5, max_order)
plt.ylim(-0.5, max_order)
adjust_spines(ax, ['left', 'bottom'])
ax.set_axisbelow(True)
plt.xlabel('$i_1$', fontsize=13)
plt.ylabel('$i_2$', fontsize=13)
cbar = plt.colorbar(extend='neither', spacing='proportional',
                orientation='vertical', shrink=0.8, format="%.0f")
cbar.ax.tick_params(labelsize=13)
plt.savefig('Pseudo_3.png',   dpi=300, bbox_inches='tight')

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.scatter(Pol.quadraturePoints[:,0], Pol.quadraturePoints[:,1] , marker='o', s=2, color='tomato')
adjust_spines(ax, ['left', 'bottom'])
plt.xlabel('$\zeta_1$', fontsize=13)
plt.ylabel('$\zeta_2$', fontsize=13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.savefig('FigureF.png' , dpi=200, bbox_inches='tight', pad_inches=0.1)

print len(np.unique(Pol.quadraturePoints, axis=0))
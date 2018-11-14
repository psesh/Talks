from equadratures import *
import numpy as np
import matplotlib.pyplot as plt


#x = Parameter(distribution='uniform', lower=-1., upper=1., order=35)
tensor = Basis('Tensor grid', [35, 35])
elements = tensor.elements

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.scatter(elements[:,0], elements[:,1], s=20, marker='o', c='navy', cmap='jet', alpha=1.0, vmin=-16.0, vmax=0.)
plt.xlim(-0.5, 35.5)
plt.ylim(-0.5, 35.5)
adjust_spines(ax, ['left', 'bottom'])
ax.set_axisbelow(True)
plt.xlabel('$i_1$', fontsize=13)
plt.ylabel('$i_2$', fontsize=13)
plt.savefig('Basis_tensor.png',   dpi=300, bbox_inches='tight')
del tensor 

tensor = Basis('Euclidean degree', [35, 35])
elements = tensor.elements

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.scatter(elements[:,0], elements[:,1], s=20, marker='o', c='navy', cmap='jet', alpha=1.0, vmin=-16.0, vmax=0.)
plt.xlim(-0.5, 35.5)
plt.ylim(-0.5, 35.5)
adjust_spines(ax, ['left', 'bottom'])
ax.set_axisbelow(True)
plt.xlabel('$i_1$', fontsize=13)
plt.ylabel('$i_2$', fontsize=13)
plt.savefig('Basis_euclidean.png',   dpi=300, bbox_inches='tight')
del tensor 

tensor = Basis('Total order', [35, 35])
elements = tensor.elements

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.scatter(elements[:,0], elements[:,1], s=20, marker='o', c='navy', cmap='jet', alpha=1.0, vmin=-16.0, vmax=0.)
plt.xlim(-0.5, 35.5)
plt.ylim(-0.5, 35.5)
adjust_spines(ax, ['left', 'bottom'])
ax.set_axisbelow(True)
plt.xlabel('$i_1$', fontsize=13)
plt.ylabel('$i_2$', fontsize=13)
plt.savefig('Basis_total_order.png',   dpi=300, bbox_inches='tight')
del tensor 

tensor = Basis('Hyperbolic basis', [35, 35], q=0.3)
elements = tensor.elements

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.scatter(elements[:,0], elements[:,1], s=20, marker='o', c='navy', cmap='jet', alpha=1.0, vmin=-16.0, vmax=0.)
plt.xlim(-0.5, 35.5)
plt.ylim(-0.5, 35.5)
adjust_spines(ax, ['left', 'bottom'])
ax.set_axisbelow(True)
plt.xlabel('$i_1$', fontsize=13)
plt.ylabel('$i_2$', fontsize=13)
plt.savefig('Basis_hyperbolic.png',   dpi=300, bbox_inches='tight')
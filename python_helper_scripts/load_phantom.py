import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from datetime import datetime
import scipy.io
import meshio




mesh   = scipy.io.loadmat('/pvfs2/Derick/EIT/Mine/data/100_samples__max_Inclusions_3__2023-07-29-14-39-04/mesh_pet.mat')
mesh
bound  = scipy.io.loadmat('/pvfs2/Derick/EIT/Mine/data/100_samples__max_Inclusions_3__2023-07-29-14-39-04/dataset_bound.mat')
domain = scipy.io.loadmat('/pvfs2/Derick/EIT/Mine/data/100_samples__max_Inclusions_3__2023-07-29-14-39-04/dataset_domain.mat')

volt = bound['outputBoundvoltage']
curr = bound['outputBoundcurrent']
angl = bound['angl_circum'       ]

cond = domain['inputConductivity']
VOLT = domain['outputVoltage'    ]

p = mesh['p']
e = mesh['e']
t = mesh['t']


print(volt.shape)
print(curr.shape)
print(angl.shape)

# print(p.shape)
# print(e.shape)
# print(t.shape)


x, y = p[0], p[1]

elements = np.delete(e, (3, 4, 5, 6), axis=0).T
triangles= np.delete(t, (3, 4, 5, 6), axis=0).T
node_vals = cond[0,:]


print(x.shape)
print(y.shape)
print(triangles.shape)
# print(node_vals.shape)
#triang = mtri.Triangulation(x, y, triangles)

print('We are plotting')
#plt.tricontourf(triang)
#plt.tricontourf(x, y, eles , node_vals, 12)
plt.colorbar()
# print('We are showing')
# plt.show()
print('We are saving')
plt.savefig('cc3.png')


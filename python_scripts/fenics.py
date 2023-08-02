from dolfin import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import scipy.io


mesh   = scipy.io.loadmat('/pvfs2/Derick/EIT/Mine/data/100_samples__max_Inclusions_3__2023-07-29-14-39-04/mesh_pet.mat')
p = mesh['p']
e = mesh['e']
t = mesh['t']
x, y = p[0].flatten()[1:10], p[1].flatten()[1:10]

number = np.shape(x)[0]

domain_vertices = []

for i in range(number):
        #print(x[i])
        domain_vertices.append(Point(x[i], y[i]))

domain = Polygon(domain_vertices)
print(domain)

mesh = generate_mesh(domain, 0.03)
plot(mesh)

print('We are saving')
plt.savefig('cc3.png')
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

# Step 1: Load the .mat file containing the mesh data
# mat_file_path = 'path/to/your/mesh_file.mat'  # Replace with the actual path to your .mat file
# data = sio.loadmat(mat_file_path)

data   = sio.loadmat('/pvfs2/Derick/EIT/Mine/data/100_samples__max_Inclusions_3__2023-07-29-14-39-04/mesh_pet.mat')
# mesh
# bound  = scipy.io.loadmat('data/1_samples__max_Inclusions_32023-07-27-19-30-27/dataset_bound.mat')
# domain = scipy.io.loadmat('data/1_samples__max_Inclusions_32023-07-27-19-30-27/dataset_domain.mat')


# Step 2: Extract the mesh data from the loaded .mat file
nodes = data['p']  # Array of node coordinates (Nx2)
elements = data['t']  # Array of element connectivity (Nx3, 3 nodes per element)
elements = np.delete(elements, (3, 4, 5, 6), axis=0).T
print(elements.shape)

# Step 3: Plot the mesh
plt.figure()
plt.triplot(nodes[:, 0], nodes[:, 1], elements-1, 'k-')  # elements-1 since MATLAB uses 1-based indexing
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Quadratic PET Mesh')
plt.gca().set_aspect('equal')
plt.savefig('cc3.png')
plt.show()

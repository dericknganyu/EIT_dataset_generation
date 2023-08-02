This directory contains

✦ src              : the code used for data generation (codes commenting for understanding to follow)
✦ dataset          : The datasets in n folders numbered parts_n. Each part is made of 5000 samples.   
✦ dataset_textures : (TBD)


In the dataset directory you would find 

➡️ mesh.mat : an object containing the mesh information

➡️ part_i (i = 1 to n) : sub-directories of 5000 samples containing

	◼️ dataset_bound.mat contains
		➤ angl_circum               (1, 339) : theta, for the point on the boundary:  
		➤ outputBoundcurrent (5000, 32, 339) : input boundary current: (maybe I should have saved it as inputBoundcurrent)
		➤ outputBoundvoltage (5000, 32, 339) : output boundary voltage: 

	◼️ dataset_domain.mat contains, 
		➤ inputConductivity   (5000, 16785) : input conductivity in the domain: 
		➤ x1, x2            each (1, 16785) : the cartesian coordinates: 
		➤ radius, theta     each (1, 16785) : the polar coordinates:  



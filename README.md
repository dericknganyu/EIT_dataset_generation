# Dataset for the paper TBD
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Python 3.85](https://img.shields.io/badge/python-3.85-blue.svg)](https://www.python.org/downloads/release/python-385/)
[![Matlab R2022a](https://img.shields.io/badge/Matlab-R2022a-orange.svg)](https://www.python.org/downloads/release/python-3100/)
[![arXiv](https://img.shields.io/badge/arXiv-xxxx.xxxxx-b31b1b.svg)](TBD)

✦ src              : the code used for data generation (codes commenting for understanding to follow)

✦ dataset          : The datasets in n folders numbered parts_n. Each part is made of 5000 samples. 

✦ dataset_textures : (TBD)


## Dataset Description

In the dataset and dataset_textures directory you would find: 

➡️ mesh.mat : an object containing the mesh information

➡️ (200) sub-directories of 100 samples each containing:

	◼️ dataset_bound.mat contains:
	
	
		➤ angl_circum              (1, 339) : theta, for the point on the boundary:  
		
		➤ outputBoundcurrent (100, 32, 339) : input boundary current: (maybe I should have saved it as inputBoundcurrent)
		
		➤ outputBoundvoltage (100, 32, 339) : output boundary voltage: 
	
	◼️ dataset_domain.mat contains: 
	
		➤ inputConductivity   (100, 16785) : input conductivity in the domain: 
		
		➤ x1, x2           each (1, 16785) : the cartesian coordinates: 
		
		➤ radius, theta    each (1, 16785) : the polar coordinates:  
	
In the src directory: 

➡️ main_v2.m        : main script used for generating the dataset

➡️ fem_eit_fwd_v2.m : finite_element solver function

➡️ gen_ellipses     : conductivity generation function

## Reference
```bibtex
@article{
}
```

## Acknowledgments




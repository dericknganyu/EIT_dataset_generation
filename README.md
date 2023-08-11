# Dataset for the paper TBD
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Python 3.85](https://img.shields.io/badge/python-3.85-blue.svg)](https://www.python.org/downloads/release/python-385/)
[![Matlab R2022a](https://img.shields.io/badge/Matlab-R2022a-orange.svg)](https://www.python.org/downloads/release/python-3100/)
[![arXiv](https://img.shields.io/badge/arXiv-xxxx.xxxxx-b31b1b.svg)](TBD)

# Runnung code
Launch matlab in parent directory run the script
```bash
main
```
# Description

In the parent directory: 

➡️ main_v2.m        : main script used for generating the dataset

➡️ fem_eit_fwd_v2.m : finite_element solver function

➡️ gen_conductivity : conductivity generation function

Sub-directories

✦ dataset               : The datasets with inclusions of uniform conductivity in folders. 

✦ dataset_textures      : The datasets with inclusions of variable conductivity in folders.

✦ python_helper_scripts : more useful python scripts and jupyter (ipython) notebooks

✦ matlab_helper_scripts : more useful matlab scripts

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
	
	◼️ inclusions.mat contains: 
	
	◼️ dataset_domain.mat contains: 
	
	◼️ 100_conductivities_samples.png: for a visualisation of samples 
	
# Example of samples with constant conductivity inclusions
![alt text](https://drive.google.com/uc?export=view&id=119M4oo_ycrYwpDIdxgvwSXwvoTc7ksTM)
[![alt text](https://github.com/dericknganyu/EIT_dataset_generation/blob/main/dataset/100_conductivities_samples.png)]

# Example of samples with variable conductivity (textured) inclusions
![alt text](https://drive.google.com/uc?export=view&id=1mtgaT4YWAA-DEyQ5AVLuKwMtDXrmtrNH)
## Reference
```bibtex
@article{
}
```

## Acknowledgmentss




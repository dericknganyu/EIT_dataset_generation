# Dataset for the paper TBD
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Matlab R2022a](https://img.shields.io/badge/Matlab-R2022a-orange.svg)](https://www.python.org/downloads/release/python-3100/)
[![Python 3.85](https://img.shields.io/badge/python-3.85-blue.svg)](https://www.python.org/downloads/release/python-385/)
[![arXiv](https://img.shields.io/badge/arXiv-2310.18636-b31b1b.svg)](https://arxiv.org/pdf/2310.18636.pdf)

## Runnung code
To generate and save ```batchSize``` samples (with ```Nmax-Nmin+1``` currents), ```numRuns``` times, 
so that each sample has not more than ```max_numInc``` inclusions, with or without textures, 
depending on the values of ```textures``` (either ```true```,```false```. or ```constant```),  navigate to the directory containing the files

➡️ In terminal, run the command
```bash
matlab -nodisplay -r "run numRuns texture max_numInc batchSize Nmax Nmin"
```

➡️ In matlab command window
```bash
run(numRuns, texture, max_numInc, batchSize, Nmax, Nmin)
```
Example: To generate the case of constant conductivity run one of the below two commands depending if in terminal or matlab command window
```bash
matlab -nodisplay -r "run 1 constant"
run(1, 'constant')
```

## FIles
In the parent directory: 

➡️ main_v2.m        : main script used for generating the dataset

➡️ fem_eit_fwd_v2.m : finite_element solver function

➡️ gen_conductivity : conductivity generation function

Sub-directories

✦ dataset               : The datasets with inclusions of uniform conductivity in folders. 

✦ dataset_textures      : The datasets with inclusions of variable conductivity in folders.

✦ dataset_constant      : The datasets with constant conductivity in folder.

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
	
### Example of samples with constant conductivity inclusions
![constant_inclusions](https://drive.google.com/uc?export=view&id=119M4oo_ycrYwpDIdxgvwSXwvoTc7ksTM)

### Example of samples with variable conductivity (textured) inclusions
![textured_inclusions](https://drive.google.com/uc?export=view&id=1mtgaT4YWAA-DEyQ5AVLuKwMtDXrmtrNH)

## Reference
```bibtex
@article{
}
```
## Acknowledgmentss





# Trajectory Planning Under Environmental Uncertainty With Finite-Sample Safety Guarantees
This repository contains the source code used for the paper titled "Trajectory Planning Under Environmental Uncertainty With Finite-Sample Safety Guarantees" published in Automatica. DOI: [10.1016/j.automatica.2021.109754](https://doi.org/10.1016/j.automatica.2021.109754)

## Dependencies
The following dependencies need to be installed/configured and must be on the MATLAB path:
- [YALMIP](https://yalmip.github.io/) (configured with an appropriate solver like [CPLEX]( https://www.ibm.com/analytics/cplex-optimizer))
- [matlab2tikz](https://github.com/matlab2tikz/matlab2tikz)
- [MagInset](https://www.mathworks.com/matlabcentral/fileexchange/49055-maginset)

## Generating plots
In order to reproduce the simulations and plots of the paper, navigate inside the inside the corresponding chapter folder and run the `generatePlots` MATLAB function. The function will run all the necessary simulations and produce all the plots of the paper. The plots will also be saved under a `plots` folder in TikZ format, which can then be readily included in a LaTeX document.

## Citing this work
Please cite the original paper when using any part of this code. BibTeX citation data:
```
@article{Lefkopoulos2021,
	author	  	= "V.~Lefkopoulos and M.~Kamgarpour",
	title	  	= "Trajectory Planning Under Environmental Uncertainty With Finite-Sample Safety Guarantees",
	journal     = "Automatica", 
	year	  	= "2021",
	volume      = "131",
	pages       = "109754",
}
```

# XBmodel_linearisation
MATLAB code to reproduce figures in the paper "Uncovering cross-bridge properties that underlie the cardiac active complex modulus using model linearisation techniques"

# Code details
Developed in MATLAB R2021b. Tested on R2017b and runs fine so will probably run on older versions of MATLAB too.  

The repository includes two scripts. The first (plot_figures.m) can be run to reproduce the figures in the results section of the paper. This calls several functions (listed in the header) as well as the data files (Saeki_data.mat and MLcolours.mat). The second script (numerical_complex_modulus.m) demonstrates the function of the complete ODE model by running through a numerical simulation of sinusodial length pertubations.

The functions XBmodel.m and XBmodel_linear.m represent the complete (ODE) implementation and the linearised implementation of the cross-bridge model presented in the paper, respectively.

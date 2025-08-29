# Super-resolved-anomalous-diffusion-code
MATLAB Code from article https://arxiv.org/pdf/2410.18133.
This code performs the fitting of the joint distribution of estimated anomalous exponent alpha and generalized diffusion coefficient D based on the assumption that the motion of particles is a fractional Brownian motion with randomly distributed physical parameters. 

For the true underlying distribution of the physical pairs {alpha_k,D_k}, two versions are possible:
1- mixture of discrete (point-like) distributions with respective weights 'c_k', that is: 

$$p(\alpha,D)=\sum\limits_{k=1}^n c_k\delta(\alpha-\alpha_k)\delta(D-D_k)$$

2- mixture of continuous distributions with respective weights 'c_k': 

$$p(\alpha,D)=\sum\limits_{k=1}^n c_k p_k(\alpha,D),$$ 

where $p_k(\alpha,D)$ is the joint PDF of the k-th component. Two user-selected forms are possible for $p_k(\alpha,D)$, in both cases $\alpha$ is assumed to have a Beta distribution over $[0,2]$, but $D$ can have either a lognormal distribution or a Rice distribution. 


The folder containing all the functions must be downloaded and added to the Matlab path. No other functions or files are required. This code has been tested on MATLAB 2021b version.

Two examples are given:


1- The file 'Example_fit_joint_PDF_alpha_D_discrete_distribution.m' is a working example for the discrete case based on simulated data with three components. 
The function 'fit_joint_PDF_alpha_D_from_trajectory.m' is doing all the work and calls all other functions.

2- The file 'Example_fit_joint_PDF_alpha_D_continuous_distribution.m' is a working example for the continuous case based on simulated data with two components. 
The function 'fit_joint_PDF_alpha_D_from_trajectory.m' is doing all the work and calls all other functions.


Required inputs are an array N x M of M trajectories of N steps, the trajectories must be given in physical units (e.g. micrometer). The frame duration dt (sec) must be given in the script.
The number of components 'n_peak' must be given.

The code gives an output with the fitted  joint distribution of estimated alpha vs ln(2D) and produces from that another figure of the empirically estimated alpha vs D as well as a histogram of estimated alpha and  another one for estimated D, all overlayed with the respective distributions of each component.

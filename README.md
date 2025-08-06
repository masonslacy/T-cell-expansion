# T-cell-expansion
Code used to simulate T cell expansion in activating micro-rod scaffolds

All code used to generate figures in "A stochastic model of T cell expansion in micro-rod scaffolds and its continuum limit: Importance of IL-2 loading and scaffold homogeneity", Lacy _et al._, is in the main script. Other required functions are in the functions file. 

Run the first section of the main script to download the required packages, and run the next section to define all default variables and initialise. Sections after this are labelled by their figure number, and running these sections will generate figures similar to those in the paper. The sections which contain code for Figs 4 and 5 also contain code to generate animations which show T cell expansion or T cell trajectories in micro-rod scaffolds.

Variables such as _T_ (which defines the total simulation time in minutes) and _ntotal_ (which is the number of simulations to average the ABM result) are set to a low value such that running the code does not require much computation time, however these values do not match those in the paper. See figures and captions in the paper to determine the correct values for _T_ and _ntotal_ if you wish to compared results to those in the paper. Other variables may be altered to generate results under different parameter regimes.

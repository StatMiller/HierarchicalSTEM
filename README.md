# HierarchicalSTEM
Code for MCMC, simulation, and figures for spatial Bayesian hierarchical model accounting for measurement error in Scanning Transmission Electron Microscope (STEM) Images as described in *Miller et al. (2019)*.

## File Descriptions
* AtomInfo.RData - PMN STEM Image data
* DataAnalysis.R - Runs MCMC on real data using 3 models (Spatial Bayesian Hierarchical Measurement Error Model, Spatial Linear Regression, and Simple Linear Regression)
* MCMC.R - MCMC for Spatial Bayesian Hierarchical Measurement Error Model
* SimStudy.R - Runs simulation study, breaking up simulations among multiple CPU's
* boxes_ancillary_functions.R - Ancillary functions used throughout other scripts, including atom column finder and MCMC's for standard models
* Figures/
  * COM_fig.R - creates diagram of A-sites, B-sites and unweighted and weighted B-site means (Figure 2 in Paper)
  * Density_Overlay.R - Creates overlaid posterior densities (Figure 3)
  * Posterior_Ellipses.R - Creates 95% posterior ellipses on subset of image (Figure 4)
  * Rect_Figures.R - Draws rectangles around atom columns and plots estimated atom column locations (Bottom Right of Figure 1)

# JAGS code

This folder contains three main reports:

## intro.Rmd

This is a very basic report that simply introduces the simplest, standard implementation of the Hierarchical Bayesian EZ-DDM. It presents **two short simulations** with fixed trial size and participant size. The main difference is that in one case we work with simulated trial-level data, while in the other we extract summary statistics directly from the proxy model's sampling distributions.

The Rmd file and its knitted HTML file can be found on the root `/JAGS/` folder. Cache results from this simulation are stored in .RData files to `/JAGS/demos/simulation-study/simplest_model`

## sample_simStudy.Rmd

This Rmd report (and its knitted HTML) can be found under `/JAGS/demos/simulation-study/design_matrix`. We present the JAGS implementation of Anne's simulation study (as shown in `/working simulation + output/model2_v2_ncp.R`).

We introduce a generalizable JAGS model that uses inner product to avoid the need to hard code the dimensions of the design matrices fed to the model. 

We showcase the flexibility of this implementation by conducting two simulation study (each with a different design) and using the same JAGS module to retrieve posterior estimates. In both cases, recovery is pretty neat.

## demo_realData.Rmd

The Rmd file (and knitted html) for this applied exercise can be found under `/JAGS/demos/applications`. It presents a STAN repplication of the applied example developed by Anne.



# A stan implemetation for the bayesian hierarchical EZ diffusion model 
*this is a work in progress, files are updated often -- for the newest updates please contact Anne Giacobello under anne.giacobello@uni-potsdam.de*

## Department of Psychology, University of Potsdam
## Authors
Anne Giacobello, Julia Haaf

## Description
This repository contains the content and code of this project. The [EZ diffusion model](https://link.springer.com/article/10.3758/BF03194023) a simplified version of the diffusion model used in cognitive psychology to analyze decision-making processes.
It simplifies the full diffusion model by focusing on three key parameters: drift rate (the speed of evidence accumulation), boundary separation (decision threshold), and non-decision time (encoding and motor response times). 
The EZ diffusion model is frequentist, and not hierarchical. Recently, there has been made an effort to [describe a probabilistic version of the EZ-diffusion model that can serve as a proxy model to the drift diffusion model](https://link.springer.com/article/10.3758/s13423-025-02729-y).
Using [JAGS](https://mcmc-jags.sourceforge.io/), they concluded that the recovery of some parameters was biased, however, generally, the recovery of (regression) parameters was good. 
We aim to put these recoveries to the test by implementing the model in [RStan](https://mc-stan.org/rstan/), once using the same priors as in the JAGS implementation and once using different priors. We also want to implement condition effects and put an emphasis on individual differences. 


## Support
For any questions, please contact Anne Giacobello (E-Mail: anne.giacobello@uni-potsdam.de). 

## Roadmap

## Project status
This project is currently under development. 

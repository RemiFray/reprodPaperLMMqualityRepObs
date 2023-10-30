# Incorporating a measure of quality and accounting for repeated observations in capture-recapture misidentification models

This repository contains all the material needed to reproduce the results from the paper.
The instructions are given below.


## Data

First generate the simulation data. In the folder "data" are 6 scripts to generate data according to models $M_{t, \alpha_n}$ or $M_{\lambda, \alpha}$ for the different numbers of occasion.
Simply run the scripts in any order.

The otter data from the study *Non-Invasive Genetic Mark-Recapture as a Means to Study Population Sizes and Marking Behaviour of the Elusive Eurasian Otter (Lutra lutra)* ([paper DOI](https://doi.org/10.1371/journal.pone.0125684)), available in a pdf ([data DOI](https://doi.org/10.1371/journal.pone.0125684.s002)) and are also in the otter.csv in the data folder.


## Run the models

The scripts used to run the models are in the folder "**mainScripts**". There are three sub-folders for:

* the simulation study with no repeated observation (simulations by the model $M_{t, \alpha_n}$), 
* the simulation study with repeated observation (simulations by the model $M_{\lambda, \alpha}$), 
* the otter study.

Most of the time, there is one script for one model and one number of occasion. For efficiency reasons, we ran each script divided in three (cheap parallelisation) but we regrouped them in here for clarity. The ones that were not regrouped are because the MCMC was not ran on the same number of iteration.  
In particular, for the model $M_{t, \alpha_n}$ the number of iteration was greatly reduced to gain time and because the reduced number was sufficient. The reason why the number of iteration may seem after a random number of simulations being analysed is that sometimes an unknown error would occur and the script stop running on the cluster it was on. We have no idea why and what the errors were. **Simply restarting the scripts from where it crashed was enough.** I believe it due to a bad interaction between NIMBLE and the cluster but if you experience such a thing, try restarting the script.

The scripts will save the markov chains in the result folder.


## Extract the results

There are script to extract summaries of the MCMC available in the folder "plotResults". There is one script per model per study. At the end of one of these scripts per study, there is also a few lines of code to merge all summary tables from the different models of the study.

There is also one script per simulation study to generate the plots in the paper.


## Functions

This folder contains all the functions, models, distributions and samplers needed to run the models. Nothing need to be touched in there.


## Results

The summary tables of results from both simulation studies and the otter study are directly available in the result folder.


# Comparison experiments between nonstationary NNGPs and state-of-the-art methods

The competitors are: 
* Nonstationary NNGP
* INLA
* Local GPs using the R package liGP
* Deep GPs using the R package deepGP

There are two rounds of experiments. 
* experiments on data sets with **nonstationary range**
* experiments on data sets with **heteroskedastic noise**
The new expermients are in the *_r3* folders.

In each of the two results folders, the *.sh* and *.R* scripts used to run the experiments on a cluster can be found. Also, a script used to aggregate the results is present. 
There also is a folder with **pdf files** corresponding to the **MCMC runs** of nonstationary NNGPs and deepGPs. 

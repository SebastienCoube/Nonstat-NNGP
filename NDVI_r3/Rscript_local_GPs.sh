#!/bin/bash -i 
 
 module load gcc/10.2.0
 module load R/4.4.0
 module load gdal/3.1.3 
 module load geos/3.11.1
 module load sqlite/3.33.0
 module load proj/7.1.1 
 

Rscript local_GPs.R

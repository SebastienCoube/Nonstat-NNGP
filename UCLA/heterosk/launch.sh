## #!/bin/bash#!/bin/bash
#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -pe shared 12
#$ -l h_rt=80:00:00,h_data=1G
#$ -M $USER@mail
#$ -m bea
 
 
 
 module load gcc/10.2.0
 module load R/4.4.0
 module load gdal/3.1.3 
 module load geos/3.11.1
 module load sqlite/3.33.0
 module load proj/7.1.1 
 
 qsub Rscript.sh 

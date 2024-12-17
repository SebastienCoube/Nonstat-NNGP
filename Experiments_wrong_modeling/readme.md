# Synthetic experiments

Presented in the annex of the paper. 

If you want to replicate those experiments, get a SLURM cluster. 
Then copy the two folders on your cluster, cd you way into one of the two folders, and use ./Run_several.sh.
This script will replicate experiment.R (into experiment1.R, experiment2.R, etc) and will launch each of them using Run.sh. 
You may want to change Run_several.sh and Run.sh to adapt the characteristics of your cluster. 
The results will appear in the foler. 
Do the same in the other folder. 
Then, rsync the results to your laptop and run Analysis_with_lm.R in each folder to get the figures of the paper (if you don't want to wait for the folders to copy, you can also run Analysis_with_lm.R on your cluster, and just pull the generated pdf files to our laptop)

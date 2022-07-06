#!/bin/bash
for i in {1..1920}; do if test ! -f "res${i}complete.RDS" ; then rm "res${i}started.RDS"  ; fi; done
for name in experiment{1..16}.R; do cp experiment.R $name; done
for name in experiment{1..16}.R; do sbatch Run.sh $name; sleep 5; done

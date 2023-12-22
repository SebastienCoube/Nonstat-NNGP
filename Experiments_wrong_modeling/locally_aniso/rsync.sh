#!/bin/bash
tatato="rsync scoube@atlas-edr.sw.ehu.es://../../scratch/scoube/aniso/res"
tatata="complete.RDS res/"
for i in `seq 1 270`;
do
echo "$tatato$i$tatata"
eval "$tatato$i$tatata"
done

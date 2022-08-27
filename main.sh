#!/bin/bash
set -e
set -u
for i in $(ls -d */)
do
cd $PWD/${i%/}
echo 'snakemake -j 25 -k --latency-wait 180 --rerun-incomplete --cluster-config ../cluster.yaml --cluster "qsub -pe OpenMP {threads} -cwd -V -l mem_free={cluster.mem} -l h_rt={cluster.time}" --jobs 25 --configfile params.yaml -s ../Snakefile' | qsub -pe OpenMP 8 -l mem_free=8G,h_rt=48:00:00 -cwd -V
cd -
done

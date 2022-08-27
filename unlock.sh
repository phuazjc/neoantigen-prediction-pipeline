#!/bin/bash
set -e
set -u
for i in $(ls -d */)
do
cd $PWD/${i%/}
snakemake -s ../tmp --unlock
cd -
done

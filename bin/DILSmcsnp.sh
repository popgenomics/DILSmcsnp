#!/usr/bin/bash
## launch DILS for 2 populations
## the provided argument is for --configfile, expecting the yaml file
module unload snakemake
module load pypy/2.7-5.10.0
#module load snakemake/5.3.0
module load snakemake/6.5.0
module load r/3.6.3
module load python/2.7
module load java-jdk/8.0.112
binpath='/home/croux/Programmes/DILSmcsnp/bin'
snakemake --snakefile ${binpath}/Snakefile_gw -p --cores 140 -j 140 --configfile ${1} --cluster-config ${binpath}/cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu}" --latency-wait 10


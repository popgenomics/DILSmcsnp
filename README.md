# DILSmcsnp  
to adapt:  
## in DILSmcsnp.sh  
line 11  
binpath='' # full link to the bin directory of the DILSmcsnp depository  
  
##in Snakefile_gw  
lines 31 and 32  
binpath='' # full path to the bin directory of the DILSmcsnp depository  
projectpath = '' # full path to the directory where all DILS analysis will be performed, i.e, where yaml files have to be placed  
  
## in yaml file  
infile: # full path to the fasta file, placed within the projectpath  
config_yaml: # full path to the current yaml file

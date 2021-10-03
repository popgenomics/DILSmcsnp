# DILSmcsnp  
to adapt to your project, the best is to set a single directory where *all* analysis (different projects for different species, etc ...) will be conducted.   
Let's name this directory *projectpath* 
The path to such *projectpath* has to be specified in the Snakefile (named here *Snakefile_gw*).  
   
## in DILSmcsnp.sh  
line 11  
*binpath=*'' # full link to the bin directory of the DILSmcsnp depository  
  
## in Snakefile_gw  
lines 31 and 32  
*binpath=*'' # full path to the bin directory of the DILSmcsnp depository  
*projectpath=* '' # full path to the directory where all DILS analysis will be performed, i.e, where yaml files have to be placed  
  
## in yaml file  
*infile:* # full path to the fasta file, placed within the projectpath  
*config_yaml:* # full path to the current yaml file

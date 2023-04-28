
``
### GREGOR

```bash


vi doGREGOR.conf
###############################################################################
INDEX_SNP_FILE = ./SNP_list
BED_FILE_INDEX = ./BED
REF_DIR = ./REF
R2THRESHOLD = 0.7 ## must be greater than 0.7
LDWINDOWSIZE = 1000000 ## must be less than 1MB; these two values define LD buddies
OUT_DIR = ./SNPlist
MIN_NEIGHBOR_NUM = 500 ## define the size of neighborhood
BEDFILE_IS_SORTED = true  ## false, if the bed files are not sorted
POPULATION = EUR  ## define the population, you can specify EUR, AFR, AMR or ASN
TOPNBEDFILES = 2
JOBNUMBER = 10
###############################################################################
#BATCHTYPE = mosix ##  submit jobs on MOSIX
#BATCHOPTS = -E/tmp -i -m2000 -j10,11,12,13,14,15,16,17,18,19,120,122,123,124,125 sh -c
###############################################################################
#BATCHTYPE = slurm   ##  submit jobs on SLURM
#BATCHOPTS = --partition=main --time=0:30:0
###############################################################################
BATCHTYPE = local ##  run jobs on local machine




--------------
vi doGREGOR.sh
--------------
#!/bin/sh
SNPlist=$1
perl GREGOR.pl --conf doGREGOR.conf
---------------------------------------------

sbatch doGREGOR.sh SNPlist





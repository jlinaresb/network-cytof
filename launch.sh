#!/bin/bash


#SBATCH -p thinnodes
#SBATCH --qos default
#SBATCH -t 10:00:00


module load cesga/2018 gcc/6.4.0 unixodbc/2.3.7 mariadb-connector-c/3.0.8 R/4.0.2
Rscript ~/git/network-cytof/parse_drugBank.r
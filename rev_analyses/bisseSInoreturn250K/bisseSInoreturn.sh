#!/bin/bash -l
#PBS -l pmem=500mb
#PBS -l nodes=1:ppn=1
#PBS -m abe
#PBS -M rzenil@umn.edu
#PBS -o bissesinoret.log
#PBS -e bissesinoreterrorscluster.log
#PBS -l walltime=48:00:00
module load gcc

/home/eeg/shared/revbayes/projects/cmake/rb  /home/eeg/shared/twostatesse/bissesinoret/bisseSInoreturn.Rev


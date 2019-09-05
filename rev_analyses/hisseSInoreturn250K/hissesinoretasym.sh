#!/bin/bash -l
#PBS -l pmem=500mb
#PBS -l nodes=1:ppn=1
#PBS -m abe
#PBS -M rzenil@umn.edu
#PBS -o hisseoutputas.log
#PBS -e hisseerrorsclusteras.log
#PBS -l walltime=64:00:00
module load gcc

/home/eeg/shared/revbayes/projects/cmake/rb /home/eeg/shared/twostatesse/hissesinoret/hissesinoretasym.Rev

#!/bin/bash -l
#PBS -l pmem=600mb
#PBS -l nodes=1:ppn=1
#PBS -m abe
#PBS -M rzenil@umn.edu
#PBS -o muhissenodipallasym.log
#PBS -e muhissenodipallasymerrorscluster.log
#PBS -l walltime=64:00:00
module load gcc

/home/eeg/shared/revbayes/projects/cmake/rb  /home/eeg/shared/twostatesse/muhissenodip/muhissenodipallasym100K.Rev
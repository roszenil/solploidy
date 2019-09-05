#!/bin/bash -l
#PBS -l pmem=500mb
#PBS -l nodes=1:ppn=1
#PBS -m abe
#PBS -M rzenil@umn.edu
#PBS -o bissenodipoutput.log
#PBS -e bissenodiperrorscluster.log
#PBS -l walltime=48:00:00
module load gcc

/home/eeg/shared/revbayes/projects/cmake/rb  /home/eeg/shared/twostatesse/bissenodip/bisseDPnodiploidization.Rev

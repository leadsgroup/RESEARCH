#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
mpiexec -n 1 /opt/python/3.7/bin/python  HC_2000ft_Mission.py
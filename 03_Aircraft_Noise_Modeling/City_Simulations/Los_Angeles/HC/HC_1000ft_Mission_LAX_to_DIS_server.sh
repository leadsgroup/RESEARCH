#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
mpiexec -n 1 /opt/python/3.7/bin/python  HC_1000ft_Mission_LAX_to_DIS_Over_Road.py

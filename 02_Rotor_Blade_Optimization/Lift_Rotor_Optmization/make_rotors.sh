#!/bin/bash
#SBATCH -J Lift_Rotor_Optimization
#SBATCH -o Rotor_Design.out
#SBATCH -e Rotor_Design.err
#SBATCH -p gpu
#SBATCH -G 2
#SBATCH -t 1-00:00 

module load python/3.9.0
module load py-numpy/1.20.3_py39
module load cuda/11.2.0
module load cudnn/8.1.1.33
module load py-scikit-learn/1.0.2_py39
module load viz
module load py-matplotlib/3.4.2_py39


python3 Lift_Rotor_Optimization.py

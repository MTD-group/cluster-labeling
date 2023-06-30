#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Abhinav Roy

Lab: MTD + McCue Research Group 
Department: Materials Science and Engineering
University: Northwestern University

script execution command:
python3 percolation.py system_size seed 
example for the case when passing command line arguments: python3 percolation.py 5 140 
"""
import fcc
import numpy as np
import analysis
import sys
import os
from time import perf_counter
from multiprocessing import Pool
# Timing the code execution
t_start = perf_counter()

# System variables as command line arguments 
sys_size = int(sys.argv[1])
seed = int(sys.argv[2])
# Defining other system variables
itr = 0
start_comp = 0.15
end_comp = 0.35
final_comp = 0
Nx = Ny = Nz = sys_size
num_atoms = 4*Nx*Ny*Nz
coord = np.zeros((num_atoms,3),dtype=float)
atom_index = []
num_comp = 6
perc_comp_list = []
# Loop for determining the percolation threshold
filename = f"realization_{seed}.txt"
f = open(os.path.join("./percolation_log", filename),"w")
f.write(f"* Realization: {seed}\n")
while (itr<10):
    itr += 1
    if itr == 1:
         composition = np.linspace(start_comp,end_comp,5,endpoint=True,dtype=float)
    else:
         composition = np.concatenate((start_comp,end_comp))
    f.write(f"Iteration: {itr}\n")
    f.write(f"Starting compositions: {composition}\n")
    comp_list = []
    pool = Pool()
    for comp in composition:
        atom_index, coord, final_comp = pool.apply(fcc.labeling, args = (Nx, Ny, Nz, num_atoms, comp, seed))
        Rx, Ry, Rz = analysis.threshold(atom_index, coord, num_atoms, Nx, Ny, Nz)
        if (Rx == 1 and Ry == 1 and Rz == 1):
                comp_list.append(final_comp)
    pool.close()
    pool.join()
    perc_comp = min(comp_list)
    f.write(f"Percolating compositions: {comp_list}\n")
    f.write(f"The system is percolating at p = {perc_comp:.4f}\n")
    perc_comp_list.append(perc_comp)
    start_comp = np.linspace(perc_comp-0.03,perc_comp,int(num_comp/2),endpoint=True,dtype=float)
    end_comp = np.linspace(perc_comp,perc_comp+0.03,int(num_comp/2),endpoint=False,dtype=float)[1:]
    f.write("\n")
f.write(f"Percolation composition list for 10 iterations: {perc_comp_list}\n")
f.write(f"! Final percolation threshold: {min(perc_comp_list):.4f}\n")
t_stop = perf_counter()
f.write(f"Time elapsed: {t_stop-t_start}")
f.close()


                                                           



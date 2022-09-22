#!/usr/bin/env python
import os
import sys
import re
from os.path import exists

with open("input.txt") as file:
    for line in file:
        line = line.rstrip()
        pdbs = line.split(',')
        for pdb in pdbs:
            path = "output_PDB/" + pdb + "_renum.pdb"
            if exists(path):
                print(pdb,"already done")
                continue
            else:
                print("running",pdb)
                os.system("python3 PDBrenum.py -rfla " + pdb + " -offz -PDB")

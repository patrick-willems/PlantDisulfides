#!/usr/bin/env python3
"""Iterate AlphaFold PDB files in folder and run DSSP, PropKa3, perform S-S distance calculations and parse corresponding PAE from the AF predictions."""

__author__  = "Patrick Willems"
__email__   = "willems533@gmail.com"           

import math
import os
import sys
import numpy
import glob
import re
import os.path
from Bio.PDB import *

pdb = sys.argv[1]
os.system('alias msms=\'/home/ubuntu/data/ss_pred/msms/msms.x86_64Linux2.2.6.1\'')
parser = PDBParser()
structure = parser.get_structure("test", pdb)
model = structure[0]
surface = get_surface(model,MSMS='/home/ubuntu/data/ss_pred/msms/msms.x86_64Linux2.2.6.1')
chain = model['A']
res52 = chain[52]
rd = residue_depth(res52, surface, MSMS='/home/ubuntu/data/ss_pred/msms/msms.x86_64Linux2.2.6.1')
print(rd)

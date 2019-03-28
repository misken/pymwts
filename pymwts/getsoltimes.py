# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 16:43:18 2012

@author: mark
"""

import glob
import os
import re

path = '/home/mark/Documents/mwts04_2/outputs/out'
files = glob.glob(path + '/*.out')
pattern1_time = re.compile("Explored ([0-9]+) nodes \(([0-9]+) simplex iterations\) in ([0-9]*\.[0-9]+|[0-9]+) seconds")
#pattern1_time = re.compile("Explored ([0-9]+) nodes \(([0-9]+) simplex iterations\)")



for f in files[0:100]:
    head, tail = os.path.split(f)
#    print tail
#    print f
    
    # Open the file
    
#Explored 1059 nodes (228963 simplex iterations) in 196.45 seconds
#Optimal solution found (tolerance 2.00e-02)
#Best objective 1.508850000000e+04, best bound 1.493197674765e+04, gap 1.0374%

    outfile = open(f,'r')
    for line in outfile:
        match1_time = pattern1_time.match(line)
        
        if match1_time:
            print tail, match1_time.group(1), match1_time.group(2), match1_time.group(3)
            break
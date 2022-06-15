#!/usr/bin/python


#seperate and clean output file of SparseDOSSA for further use 

import sys
import os
import timeit
import glob
import re

import numpy as np
from collections import Counter

def read_file(fname):
    #----------------- read input -----------------------
    with open(fname) as ff:
        con = ff.readlines()
        for line in ff:
            data = line.split()
            print(data)


    ff.close()
    #----------------------------------------------------
    return con


#---------read input -----------------------
ft = 'otu_table.tsv'
con = []
con =read_file(ft)

#split inputfile into three files
for i in con:
    designation = i. split("\t")[1]
 
    if (designation.startswith( 'Feature_Lognormal_' )):
        lognorm = open("lognorm_t.tsv", "a")
        lognorm.write(i)
    if (designation.startswith( 'Feature_Outlier_' )):   
        outlier = open("outlier_T.tsv", "a")
        outlier.write(i)
    if (designation.startswith( 'Feature_spike_n_1_m_1_')):
        spike = open("spike_t.tsv", "a")
        spike.write(i)
    




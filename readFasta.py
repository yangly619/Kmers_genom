#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Time  : 04/05/2019 12:10
# @Author: Yangyang Li
# @File  : readFasta.py

#from Bio import SeqIO

import pandas as pd
import itertools as it

def read(file):
    f = open(file, "r")
    list = dict()
    accNumber = []
    species = []
    seq = []
    for line in f:
        if line.startswith('>'):
            line=line.split(" ",maxsplit=1)
            accNumber.append(line[0])
            species.append(line[1].split(","))
        elif line.startswith('A') or line.startswith('C') or line.startswith('G') or line.startswith('T'):
            seq.append(line.strip("\n"))

    list = pd.DataFrame([accNumber, species, seq], columns=["accNum","species", "sequence"])

    return list
print(read("secuencias.txt"))

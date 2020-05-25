#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Time  : 04/05/2019 12:10
# @Author: Yangyang Li
# @File  : readFasta.py
import sys
import time
from logging import Logger

from Bio import SeqIO
import kmerCalculator as kc
from pandas.core.frame import DataFrame

def read(file):
    f = open(file, "rU")
    seq = {}
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            info = line.replace('>', '')
            seq[info] = ''
        else:
            seq[info] += line.replace('\n', '').strip()
            # seq[info] += line.strip('\n')
    f.close()
    return seq


def load_seq(file):
    seq_list = list()
    with open(file, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_list.append(record.seq)
    return seq_list


def load_id(file):
    id_list = list()
    with open(file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id_list.append(record.id)
    return id_list


def load_info(file):
    info_list = list()
    with open(file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            info_list.append(record.description)
    return info_list


def freq_table(file,k,alp):
    kmer = []
    for i in range(0, 4 ** 2):
        kmer.append(kc.kmer_index2word(i, k, alp))
    freq= []
    with open(file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            freq.append(kc.compute_freq(kc.compute_kmers(record.seq,k,alp)))
    return DataFrame(freq,columns=kmer)

#print(load_id("lactobacillus.fasta"))
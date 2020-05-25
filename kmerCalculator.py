#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Time  : 04/05/2019 13:42
# @Author: Yangyang Li
# @File  : kmerCalculator.py

def kmer_index(sequence, k, alphabet):
    res = 0
    if len(sequence) < k:
        return -1
    for i in range(0, k):
        res += alphabet.index(sequence[i]) * 4 ** (k - i - 1)
    return res


def change_char(string, position, replacement_char):
    return string[:position] + replacement_char + string[position + 1:]


def kmer_index2word(kmer_hash, kmer_size, alphabet):
    quotient = kmer_hash
    word = 'A' * kmer_size
    i = 0
    while quotient >= 4:
        reminder = quotient % 4
        quotient = int(quotient / 4)
        word = change_char(word, kmer_size - i - 1, alphabet[reminder])
        i += 1

    word = change_char(word, kmer_size - i - 1, alphabet[quotient])
    return word


################### version1

def calculate_Frecuency(seq, k, alphabet, list):
    indice = kmer_index(seq, k, alphabet)
    list[indice] = list[indice] + 1
    total = sum(list)
    list[indice] = round(list[indice] / total, 2)
    return total


def compute_kmer(k, sequence, alphabet):
    seq = 0
    l = []
    lst = [0 for n in range(4 ** k)]
    for i in range(len(sequence)):
        seq = sequence[i:i + k]
        if len(seq) == k:
            calculate_Frecuency(seq, k, alphabet, lst)
    return lst


###########################version2

def compute_kmers(sequence, kmer_size, alphabet):
    frequency = [0] * 4 ** kmer_size
    i = 0
    while i < len(sequence) - (kmer_size - 1):
        n = kmer_index(sequence[i:i + kmer_size], kmer_size, alphabet)
        frequency[n] += 1
        i += 1
    return frequency


def compute_freq(kmers):
    total = sum(kmers)
    for i in range(0, len(kmers)):
        kmers[i] = round(kmers[i] / total, 4)

    return kmers

#print(compute_freq(compute_kmers("AAAAC",2,"ACGT")))
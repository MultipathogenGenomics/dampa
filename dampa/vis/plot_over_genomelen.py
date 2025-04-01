from statistics import median
import random
import pandas as pd
import os
import time
import sys
from datetime import timedelta
import glob
import re
import pandas as pd
from fontTools.cffLib.specializer import programToString
from matplotlib import pyplot as plt
import seaborn as sns
from itertools import groupby
from collections import Counter
from time import sleep as sl

def count_zero_sequences(lst):
    # Group consecutive elements
    groups = [(k, len(list(g))) for k, g in groupby(lst)]

    # Extract only the lengths of sequences of 0s
    zero_lengths = [length for value, length in groups if value == 0]

    # Count occurrences of each sequence length
    return Counter(zero_lengths)

def replace_short_zeros(lst,threshold):
    """
    function needed becuase blast will call 0 at ends of probe if there are snps (due to it being a local alignment) but capture will not be effected by this
    """
    n = len(lst)
    i = 0
    inls = list(lst)
    while i < n:
        if lst[i] == 0:
            start = i
            while i < n and lst[i] == 0:
                i += 1
            if i - start < threshold:  # If the sequence of 0s is shorter than threshold
                for j in range(start, i):
                    lst[j] = 1
        else:
            i += 1
    return lst

def load_capture_data(capture_path):
    # print(f'Loading capture data from {capture_path}...')
    with open(capture_path, 'r') as input_file:
        headers, seqs, depths = [], [], []
        for line in input_file:
            if line[0] == '>':
                if not seqs == [] and seqs[-1] == '':
                    print(header[-1])
                header = line.strip().lstrip('>')
                headers.append(header)
                seqs.append('')
                depths.append('')
            elif line[0] == '$':
                seqs[-1] += line.strip()
            elif line[0] == '#':
                depths[-1] += line.strip()
    seqs = [seq.lstrip('$') for seq in seqs]
    depths = [[int(d) for d in depth.lstrip('#').split(',')] for depth in depths]
    depths = [replace_short_zeros(x, 10) for x in depths]
    if len(set(len(c) for c in [headers, seqs, depths])) != 1:
        # logger.error(f'The number of headers, seqs, and probe depth lists do not match in {capture_path}.')
        exit(1)
    for header in headers:
        if headers.count(header) > 1:
            # logger.error(f'Header {header} appears more than once in {capture_path}.\n')
            exit(1)
    capture_data = {header: (seq, depth) for header, seq, depth in zip(headers, seqs, depths)}
    for header in capture_data.keys():
        if len(capture_data[header][0]) != len(capture_data[header][1]):
            # logger.error(f'Seq length and probe depth list length do not match for entry {header}.')
            exit(1)
    # print(f' Total targets loaded: {"{:,}".format(len(capture_data))}')
    return capture_data



def filt_lens(capdata,headers=False):
    lens =[len(capdata[x][1]) for x in capdata]
    medlen = median(lens)
    medtop = medlen + (medlen*0.1)
    medbottom = medlen - (medlen*0.1)

    filtlensd = {x:capdata[x] for x in capdata if len(capdata[x][1]) < medtop and len(capdata[x][1]) > medbottom}

    percpos = []
    filtlens = []
    hfilt = []
    for k,v in filtlensd.items():
        dep=v[1]
        filtlens.append(dep)
        pp = [(x / len(dep)) * medlen for x in range(1, len(dep) + 1)]
        percpos.append(pp)
        hfilt.append(k)
    zipped = list(zip(filtlens,percpos,hfilt))
    return zipped

def make_genome_plots(capdata,pref):
    # print("load")
    zipdata = filt_lens(capdata)
    # print("normalise")
    totno = len(zipdata)
    if totno > 1000:
        zipdata = random.sample(zipdata,800)
        totno = 800
    alphaset = float(10)/totno

    c=0
    plt.figure(figsize=(10, 2))
    for data in zipdata:
        # print(name)
        depths = data[0]
        # make_coloured_line(range(len(depths)),depths,alphaset)
        sns.lineplot(x=range(len(depths)),y=depths,alpha=alphaset,color='blue')
    plt.axhline(y=1, color='black', linestyle='--',zorder=0,alpha=0.3)
    plt.savefig(f"{pref}_genomecov.pdf")


    plt.clf()

    plt.figure(figsize=(10, 2))
    c=0
    for data in zipdata:
        # name=data[2]
        depths = data[0]
        proplen = data[1]
        try:
            # make_coloured_line(proplen, depths, alphaset)
            sns.lineplot(x=proplen,y=depths,alpha=alphaset,color='blue')
        except:
            print(len(depths))
            print(len(proplen))
        c+=1
        # if c%1000 == 0:
        #     print(f"p2{c}")
    plt.axhline(y=1, color='black', linestyle='--',zorder=0,alpha=0.3)
    plt.savefig(f"{pref}_genomecov_lennorm.pdf")

if __name__ == "__main__":
    inpt = ""
    outpref = ""

    capdata = load_capture_data(inpt)
    make_genome_plots(capdata,outpref)
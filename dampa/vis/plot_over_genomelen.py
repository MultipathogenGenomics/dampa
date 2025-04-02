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


def replace_short_zeros(lst,threshold):
    """
    function needed becuase blast will call 0 at ends of probe if there are snps (due to it being a local alignment) but capture will not be effected by this
        Args:
        lst (list): List of integers.
        threshold (int): Threshold length for replacing zeros.

    Returns:
        list: Modified list with short sequences of zeros replaced by ones.
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


def filt_lens(capdata, headers=False):
    """
    Filters capture data based on the median length of sequences.

    Args:
        capdata (dict): Dictionary containing capture data.
        headers (bool): Whether to include headers in the output.

    Returns:
        list: List of tuples containing filtered lengths, percentage positions, and headers.
    """
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

def make_genome_plots(capdata, pref):
    """
    Generates genome coverage plots from capture data.

    Args:
        capdata (dict): Dictionary containing capture data.
        pref (str): Prefix for output files.

    Returns:
        None
    """
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
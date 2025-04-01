"""
## Author: Michael Payne, michael.payne@sydney.edu.au

This script takes
"""
import argparse
import glob
import logging
import time
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def process_count_pt(countfile, maxcheck):
    countdf = np.array(countfile[1])
    propcovs = []
    # countdf["depth"] = countdf.iloc[:, columns_to_sum].sum(axis=1)
    nnocov = [countfile[1][x] for x in range(len(countfile[1])) if countfile[1][x] == 0 and countfile[0][x] == "N"]
    countnnocov = len(nnocov)
    for c in range(maxcheck + 1):
        count_over_threshold = len([x for x in countfile[1] if x == c])
        if c == 0:
            count_over_threshold = count_over_threshold - countnnocov
        total_rows = len(countfile[1]) - countnnocov
        propcov = (count_over_threshold / total_rows) * 100
        propcovs.append(propcov)
    count_over_threshold = len([x for x in countfile[1] if x > maxcheck])

    # nnocov = [countfile[0][x] for x in range(len(countfile[1])) if countfile[1][x]==0]

    if countnnocov > 5:
        ...
    total_rows = len(countdf) - countnnocov
    propcov = (count_over_threshold / total_rows) * 100
    propcovs.append(propcov)
    countdf = np.array(countfile[1])
    meancov = countdf.mean()

    return propcovs,meancov,countfile[1]


def make_propsplot(propcovs,names,maxcheck,outpath,clusters=False):
    inls = {names[x]:propcovs[x] for x in range(len(propcovs))}
    indata = pd.DataFrame(inls).T
    indata.columns = [f"{x}" for x in range(maxcheck+1)] + [f">{maxcheck}","meancov"]
    if clusters:
        indata["cluster"] = indata.index.map(clusters)
        indatamelt = indata.reset_index()
        indatamelt = indatamelt.melt(id_vars=['index', 'cluster'],var_name="Sample", value_name="Values")
        indatamelt = indatamelt[indatamelt["Sample"] != "meancov"]
        plt.figure(figsize=(8, 6))
        plt.title("Percentage of each genome covered by a given probe depth")
        plt.xlabel("Probe depth")
        plt.ylabel("Percentage of genome covered")
        plt.ylim(-5, 100)
        sns.boxplot(data=indatamelt, x="Sample", y="Values", hue='cluster', palette='tab20')
        plt.tight_layout()
        plt.savefig(f"{outpath}_perc_genome_cov_boxplots.pdf")
        plt.clf()
        plt.figure(figsize=(8, 6))
        plt.title("Percentage of each genome covered by a given probe depth")
        plt.xlabel("Probe depth")
        plt.ylabel("Percentage of genome covered")
        plt.ylim(-5, 100)
        plot = sns.stripplot(data=indatamelt, x="Sample", y="Values", jitter=0.4, palette='tab20', hue="cluster", zorder=1)

        handles, labels = plot.get_legend_handles_labels()
        category_counts = indatamelt.loc[indatamelt["Sample"]=='0']['cluster'].value_counts()
        new_labels = [f"{label} ({category_counts[label]})" if label in category_counts else label for label in labels]
        plot.legend(handles, new_labels, title='Cluster with Counts')


        plt.tight_layout()
        # plt.show()
        plt.savefig(f"{outpath}_perc_genome_cov_stripplot.pdf")
    else:
        indatamelt = indata.reset_index()
        indatamelt = indatamelt.melt(id_vars=['index'],var_name="Sample", value_name="Values")
        indatamelt = indatamelt[indatamelt["Sample"] != "meancov"]
        plt.figure(figsize=(8, 6))
        plt.title("Percentage of each genome covered by a given probe depth")
        plt.xlabel("Probe depth")
        plt.ylabel("Percentage of genome covered")
        plt.ylim(-5, 100)
        sns.stripplot(data=indatamelt, x="Sample", y="Values", jitter=0.4, hue="Sample",palette='Set2', zorder=1)
        sns.boxplot(data=indatamelt, x="Sample", y="Values",boxprops=dict(facecolor='none', linewidth=2),
                medianprops=dict(color='black', linewidth=2),zorder=2,showfliers=False)
        plt.tight_layout()
        plt.savefig(f"{outpath}_perc_genome_cov_boxplots.pdf")


def make_stats(propcovs,names,maxcheck,meancovs,outpath,probes,report0covperc):
    inls = dict({names[x]:propcovs[x] for x in range(len(propcovs))})
    [inls[names[x]].append(meancovs[x]) for x in range(len(propcovs))]

    indata = pd.DataFrame(inls).T
    indata.columns = [f"{x}" for x in range(maxcheck+1)] + [f">{maxcheck}","meancov"]
    failcovcheck = (indata['0'] > report0covperc).sum()
    indatamean = indata.mean(axis=0)
    indatamin = indata.min(axis=0)
    indatamax = indata.max(axis=0)
    indatamedian = indata.median(axis=0)
    stats=pd.concat([indatamean,indatamin,indatamax,indatamedian],axis=1)
    stats.columns=["mean","min","max","median"]
    new_row = pd.DataFrame([{"mean": probes}],index=["total_probes"])
    stats = pd.concat([stats, new_row])
    new_row = pd.DataFrame([{"mean": failcovcheck}],index=["failed_genomes"])
    stats = pd.concat([stats, new_row])
    indata.to_csv(f"{outpath}_pergenome_stats_of_perc_coverage.csv",float_format='%.3f')
    stats.to_csv(f"{outpath}_summary_stats_of_perc_coverage.csv",float_format='%.3f')






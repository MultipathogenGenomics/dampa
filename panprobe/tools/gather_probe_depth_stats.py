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

class RuntimeFormatter(logging.Formatter):
    def __init__(self, fmt=None, datefmt=None):
        super().__init__(fmt, datefmt)
        self.start_time = time.time()  # Record the start time

    def format(self, record):
        # Calculate elapsed time
        elapsed_time = time.time() - self.start_time
        record.runtime = f"{elapsed_time:.2f}s"  # Add runtime to the record
        return super().format(record)

def setup_logging():
    formatter = RuntimeFormatter("%(asctime)s - [Runtime: %(runtime)s] - %(levelname)s - %(message)s")
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    logger = logging.getLogger("runtime_logger")
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    return logger

def get_countfiles(input,logger):

    files = glob.glob(f"{input}/*pileupcount.txt")
    if len(files) == 0:
        logger.error("No count files found, please check your input paths.\n Inputs must end in *pileupcount.txt")
        sys.exit(1)
    else:
        logger.info(f"Found {len(files)} count files")
    return files

def process_count(countfile,maxcheck):
    countdf = pd.read_csv(countfile, sep=",")
    propcovs = []
    # sum columns for ATCGN to get depth
    columns_to_sum = [2, 3, 4, 5, 7]
    countdf["depth"] = countdf.iloc[:, columns_to_sum].sum(axis=1)

    for c in range(maxcheck+1):
        count_over_threshold = (countdf["depth"] == c).sum()
        total_rows = len(countdf)
        propcov = (count_over_threshold / total_rows) * 100
        propcovs.append(propcov)
    count_over_threshold = (countdf["depth"] > maxcheck).sum()
    total_rows = len(countdf)
    propcov = (count_over_threshold / total_rows) * 100
    propcovs.append(propcov)

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

def make_allboxplot(depthls,names,args):
    df = pd.DataFrame(depthls).transpose()
    df.columns = names

    # Melt the DataFrame to long format for seaborn
    df_melted = df.melt(var_name="Sample", value_name="Values")

    plt.figure(figsize=(len(depthls) * 0.15, 6))
    sns.stripplot(x="Sample", y="Values", data=df_melted, linewidth=0.2, size=2, color="lightblue",jitter=True)
    means = df.mean()
    for i, mean in enumerate(means):
        plt.plot(i, mean, marker="o", markersize=2, color='red', label="Mean" if i == 0 else "")
    plt.axhline(y=1, color='red', linestyle='-', linewidth=2, zorder=-1)
    #plt.figure(figsize=(len(depthls) * 0.1, 6))
    # flierprops = dict(marker='.', color='black', markersize=1,markerfacecolor='black')
    # plt.boxplot(depthls,labels=names,flierprops=flierprops, widths=0.2)
    # plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
    plt.xlabel("genome name")
    plt.xticks([])
    plt.ylabel("Probe depth")
    plt.title("Boxplot of probe depth for all samples")
    plt.tight_layout()
    plt.show()

def make_1dheatmap(depthls,names,outpath):

    indict = {names[x]:depthls[x] for x in range(len(depthls))}
    max_len = max(len(v) for v in indict.values())

    # Ensure all lists have the same length by padding with NaN
    for key, value in indict.items():
        if len(value) < max_len:
            indict[key].extend([np.nan] * (max_len - len(value)))

    df = pd.DataFrame.from_dict(indict, orient='columns')

    # Convert to DataFrame of counts
    counts_df = df.apply(lambda col: col.value_counts()).fillna(0).astype(int)
    counts_df = counts_df[counts_df.index <= 20]

    greater_than_20_row = df.apply(lambda col: (col > 20).sum())  # Counts the values greater than 20 in each column
    greater_than_20_row.name = '>20'# Set the row name as '>20'
    greater_than_20_row_df = greater_than_20_row.to_frame().T
    # Append the '>20' row to the counts_df
    counts_df = pd.concat([counts_df, greater_than_20_row_df])
    counts_df = counts_df.div(df.sum(axis=0), axis=1) * 100
    # Create a vertical 1D heatmap for each column
    plt.figure(figsize=(len(counts_df.T)*0.4, len(counts_df)*0.2))  # Set figure size (taller for more categories)

    mask = counts_df == 0
    # Define a custom colormap with white for zero counts
    cmap = sns.color_palette("viridis", as_cmap=True)
    cmap.set_under("white")  # Set color for values under the min (zero in this case)

    ax = sns.heatmap(counts_df, cmap="viridis_r", cbar=True, linewidths=0,mask=mask)

    # Make only vertical lines visible to make column effect
    for i in range(counts_df.shape[0] + 1):
        ax.axhline(i, color='white', lw=0)
    for i in range(counts_df.shape[1] + 1):
        ax.axvline(i, color='white', lw=4)

    plt.gca().invert_yaxis()

    plt.title('Percentage of each genome at each probe depth')
    plt.ylabel('Probe depth')
    plt.xlabel('Genome name')
    plt.tight_layout()
    plt.savefig(f"{outpath}_sample_depth.pdf")

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



def make_meanplot(meancovs,names,ouput):
    ...



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
    indata.to_csv(f"{outpath}pergenome_cov_stats.csv",float_format='%.3f')
    # stats.loc[len(stats)] = [len(mapped), None, None, None]
    # stats.loc[len(stats)] = ["mapped", None, None, None]
    # inls["mapped"] = [len(mapped),"","","","","",""]
    # inls["total"] = [len(mapped)+len(unmapped), "", "", "","","",""]

    stats.to_csv(f"{outpath}summary_stats_of_perc_coverage.csv",float_format='%.3f')



def main():
    """
    iterate over pileupcount files in input
    get overall depth at each position
    for each genome get
        average depth
        % of genome with 0
        % of genome with >1
        % of genome with 1
    for whole dataset plot distribution of coverage averages
    could also plot boxplot for each genome with range of depths
    """
    global logger
    logger = setup_logging()
    logger.info("Starting")

    args = get_args()
    logger.info("Loading data")
    countfiles = get_countfiles(args.input,logger=logger)
    depthls = []
    names = []
    meancovs = []
    propcovs = []
    maxcheck = 4
    for countfile in countfiles:
        name=countfile.split("/")[-1].replace("pileupcount.txt","")
        stats = process_count(countfile,maxcheck) # propcovs,meancov,alldepth
        propcovs.append(stats[0])
        depthls.append(stats[2])
        names.append(name)
        meancovs.append(stats[1])
    logger.info("Making boxplots")
    # make_allboxplot(depthls,names,args)
    # make_1dheatmap(depthls,names,args)
    mapped = set()
    unmapped = set()
    make_stats(propcovs,names,maxcheck,meancovs,args,mapped,unmapped)
    make_propsplot(propcovs,names,maxcheck,args)
    make_meanplot(meancovs,names,args)

    logger.info("Finished")


def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-i", "--input", required=True, type=str, help="Folder containing input depth files produced by shiver/bin/tools/AnalysePileup.py")
    # parser.add_argument("-o", "--output", required=True, type=str,
    #                     help="output summary file")
    args = parser.parse_args()
    args.input = ""
    # args.input = "testing/"
    args.output = "testing/testout_all"
    return args

if __name__ == "__main__":

    main()



# A script for selecting unique genomes from a set of redundant genomes 
# based on pairwise comparisons and contamination and competion estimates
# to pick the best representatives
##########################################################################
#Usage:
# requires python3, networkx, operator 
#
# python selectUniqueGenomes.py pariWiseTable binstats.txt ANI COV
#
##########################################################################
import sys
import networkx as nx
from networkx.algorithms.components.connected import connected_components
import matplotlib.pyplot as plt; plt.ion()
import operator

class Bin:
    def __init__(self,name):
        self.name = name
        self.comp = -1
        self.cont = -1
        self.comparisons = list()

class Comparison:
    def __init__(self,name,ani,cov):
        self.compared_to = name
        self.cov = cov
        self.ani = ani

    def __str__(self):
        return (" ".join([self.compared_to,self.ani,self.cov]))

    def pass_threshold(self,ani_thresh,cov_thresh):
        if ani >= ani_thresh and cov >= cov_thresh:
            return True
        else:
            return False


def read_comparisons(pairwise_table,ani_thresh,cov_thresh):
    genome_list = set()
    bin_list = list()

    with open(pairwise_table) as comparisons:
        for comparison in comparisons:
            row,bin1,bin2,ani,cov = comparison.split("\t")
            temp_comparison = Comparison(bin2,ani,cov)
            if temp_comparison.pass_threshold(ani_thresh,cov_thresh):
                print(temp_comparison)





############################### MAIN #####################################################
ani = .99 #float(sys.argv[3])
cov = .75 #float(sys.argv[4])
pairwise_file = "C:\\Users\\jteva\\Git-Projects\\Pairwise_Dereplication\\pairwise_long.tsv"#sys.argv[1]
checkm_file = "C:\\Users\\jteva\\Git-Projects\\Pairwise_Dereplication\\checkm_results.tsv"#sys.argv[2]

#Read file into comparisons and make genome set
read_comparisons(pairwise_file,ani,cov)
#Select set of comparisons that meet requirements
#Create bins from comparisons (only bins that meet threshold reqs)
#create distance matrix with bins


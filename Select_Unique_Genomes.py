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

class Comparison:
    def __init__(self,a,b):
        self.comp_forward = str(a + ";" + b)
        self.comp_reverse = str(b + ";" + a)
        self.bin_a = a
        self.bin_b = b
        self.cov_forward = -1
        self.cov_reverse = -1
        self.sim_forward = -1
        self.sim_reverse = -1

    def __str__(self):
        return (self.comp_forward+" "+ str(self.sim_forward) +" "+ str(self.cov_forward) +
               "\n"+self.comp_reverse +" "+ str(self.sim_reverse) +" "+ str(self.cov_reverse))

    def fillData(self,pairwiseDataDict):
        found_forward = pairwiseDataDict[self.comp_forward]
        self.cov_forward = float(found_forward[1])
        self.sim_forward = float(found_forward[0])

        found_reverse = pairwiseDataDict[self.comp_reverse]
        self.cov_reverse = float(found_reverse[1])
        self.sim_reverse = float(found_reverse[0])

class RedundantGroup(set):
    def __init__(self,genome):
        self.genome = genome

def readData(inFile):
    genomeSet = set()
    compDict = dict()
    with open(inFile) as pairTable:
        for line in pairTable:
            if "bin1" not in line:
                splitLine = line.split("\t")
                genomeSet.add(splitLine[1])
                genomeSet.add(splitLine[2])
                compDict[';'.join([splitLine[1],splitLine[2]])] = [splitLine[3],splitLine[4]]
    genomeSet = list(genomeSet)
    return genomeSet,compDict

def readBinStats(binStatsFile):
    binStats = dict()
    with open(binStatsFile) as stats:
        for stat in stats:
            if "Marker" not in stat:
                splitStat = stat.split("\t")
                #compute bin score completeness - 5*contamination
                score = float(splitStat[11]) - (float(splitStat[12]) * 5)
                binStats[splitStat[0]]= score
    return binStats

def makeRedundantGroups(genomeList,comparisonList):
    redundantSetList = list()
    for i in genomeList:
        tempCompList = [x for x in comparisonList if i in x.comp_forward]
        tempGroup = RedundantGroup(i)
        for j in tempCompList: 
            tempGroup.add(j.bin_a)
            tempGroup.add(j.bin_b)
        redundantSetList.append(tempGroup)
    return(redundantSetList)

def makeGraph(redundantSetList): 
    G = nx.Graph()
    for genomeSet in redundantSetList:
        # each sublist is a bunch of nodes
        G.add_nodes_from(genomeSet)
        # it also imlies a number of edges:
        G.add_edges_from(makeEdges(genomeSet))
    return G

def makeEdges(set):

    it = iter(set)
    last = next(it)

    for current in it:
        yield last, current
        last = current

def plotGraph():
    pass

def selectBest(uniqGroupsList,binStats):
    representatives = list()
    for group in uniqGroupsList:
        subBinStats = dict((x, binStats[x]) for x in group if x in binStats)
        groupRep = max(subBinStats.items(), key=operator.itemgetter(1))[0]
        representatives.append(groupRep)
    return representatives

############################### MAIN #####################################################
ANI=float(sys.argv[3])
COV=float(sys.argv[4])
genomeSetList,compDict = readData(sys.argv[1])
#genome set is uniq list of bins in pairwise comparisons
#compDict is data for each comparison in pairwise comparisons

binStats = readBinStats(sys.argv[2])

#init comparisons
comparisonsList = list()
for i in range(0, len(genomeSetList), 1):
    for j in range(0, len(genomeSetList), 1):
        comparisonsList.append(Comparison(genomeSetList[i],genomeSetList[j]))

#fill comparison objects with data
for i in comparisonsList:
    i.fillData(compDict)

#subset to comparisons that have 99% ident for both forward and rev, cov of > 75% and not a comparison to itself
redundantComparisons = [x for x in comparisonsList if x.sim_forward > ANI and x.sim_reverse > ANI
                    and x.cov_forward > COV and x.cov_reverse > COV and x.comp_forward != x.comp_reverse]


#make list of RedundantGroups
redundantSetList = makeRedundantGroups(genomeSetList,redundantComparisons)
#find any sets that share entries and merge them to create a
#new set that is unique
nonUniqGroups = [x for x in redundantSetList if len(x) >0]
G = makeGraph(nonUniqGroups)
uniqGroups = connected_components(G)


#select best representative from each uniq set with dparks metric
#then add them to unique genomes
groupRepresentatives = selectBest(list(uniqGroups),binStats)

# make final list by adding originally uniq genomes and found representatives
uniqueGenomes = ([x.genome for x in redundantSetList if len(x) == 0] +  
                         groupRepresentatives)

# write list to file
outfileName = "unique_genomes_"+str(ANI)+"ani_"+str(COV)+"cov.txt"
outfile = open(outfileName,'w+')
for genome in uniqueGenomes:
    outfile.write(genome + "\n")
outfile.close()

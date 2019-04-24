import pandas as pds
import networkx as nx
import operator 
import click

class Genome:
    def __init__(self,name):
        self.name = name
        self.qual_score = -1
        self.comparison_info = dict()
    
    def __str__(self):
        return ("name: " + self.name + "\n" +
                "qual_score: " + str(self.qual_score) + "\n" +
                "comparison Count: " + str(len(self.comparison_info)) + "\n")

def trim_ext(bin_name):
    if ".fa" in bin_name:
        return bin_name.replace(".fa","")
    if ".fasta" in bin_name:
        return bin_name.replace(".fasta","")
    if ".fna" in bin_name:
        bin_name.replace(".fna","")
    else:
        return bin_name
    
def create_genome_list(pairwise_tsv,checkm_tsv):
    genome_list = list()
    with open(pairwise_tsv) as comparisons:
        for comparison in comparisons:
            if "bin1" not in comparison:
                row_num,bin1,bin2,ani,cov = comparison.strip().split("\t")
                bin1 = trim_ext(bin1)
                bin2 = trim_ext(bin2)
                if len([x for x in genome_list if x.name == bin1]) == 0:
                    new_genome = Genome(bin1)
                    new_genome.comparison_info[bin2] = [float(ani),float(cov)]
                    genome_list.append(new_genome)
                else:
                    existing_genome = [x for x in genome_list if x.name == bin1][0]
                    existing_genome.comparison_info[bin2] = [float(ani),float(cov)]
    with open(checkm_tsv) as binstats:
        for stat in binstats:
            if "Completeness" not in stat:
                split_stat = stat.strip().split("\t")
                bin_name = trim_ext(split_stat[0])
                bin_score = float(split_stat[11]) - (float(split_stat[12]) * 5)

                genome_to_update = [x for x in genome_list if x.name == bin_name][0]
                genome_to_update.qual_score = bin_score
                
    return(genome_list)

def create_matrix(genome_list,value_type):
    if value_type == "ani":
        index = 0
    if value_type == "cov":
        index = 1
        
    matrix = dict()
    for genome in genome_list:
        value_dict = dict()
        for key, value in genome.comparison_info.items():
            value_dict[key] = value[index]
        matrix[genome.name] = value_dict
    return matrix

def selectBest(uniqGroupsList,genome_list):
    representatives = list()
    for group in uniqGroupsList:
        #get genomes in the group and compare score
        genome_sub = [x for x in genome_list if x.name in group]
        representatives.append(max(genome_sub, key=lambda genome: genome.qual_score).name)
    return representatives

@click.command()
@click.option('--ani',default=.99,help="An ani value between 0 and 1 to use as the threshold for creating clusters")
@click.option('--cov',default=.75,help="An cov value between 0 and 1 to use as the threshold for creating clusters")
@click.option('--genome_comparisons',help="A tab delimited file of genome comparisons with the columns: bin1,bin2,ANI,Cov")
@click.option('--genome_stats',help="A tab delimited file of genome stats from checkM (can be produced by using checkms --tab_table option)")
def dereplicate(ani,cov,genome_comparisons,genome_stats):
    ANI = ani
    COV = cov
    PAIRWISE_TSV = genome_comparisons
    CHECKM_TSV = genome_stats

    genome_list = create_genome_list(PAIRWISE_TSV,CHECKM_TSV)

    genome_graph = nx.MultiGraph()
    genome_graph.add_nodes_from([x.name for x in genome_list])

    for genome in genome_list:
        for comparison in genome.comparison_info:
            temp = genome.comparison_info[comparison]
            ani = temp[0]
            cov = temp[1]
            if ani >= ANI and cov >= COV and genome.name != comparison:
                genome_graph.add_edge(comparison,genome.name,weight=ani)
            
    #trim edges from clusters that dont have pairwise comparison meeting threshhold
    adj_iter = genome_graph.adjacency()
    remove_list = list()
    for node in adj_iter:
        node_name = node[0]
        neighbors = node[1]
    
        for key in neighbors:
            if len(neighbors[key]) < 2:
                remove_list.append((node_name,key))

    genome_graph.remove_edges_from(remove_list)    

    #get clusters and unique genomes
    clusters = [x for x in nx.connected_components(genome_graph) if len(x) > 1]
    unique_genomes = [list(x)[0] for x in nx.connected_components(genome_graph) if len(x) == 1]

    #select representative from each cluster
    reps_from_clusters = selectBest(clusters,genome_list)
    final_bins = unique_genomes + reps_from_clusters

    outfileName = "unique_genomes_"+str(ANI)+"ani_"+str(COV)+"cov.txt"
    with open(outfileName,"w") as outfile:
        outfile.write("\n".join(final_bins))

if __name__ == '__main__':
    dereplicate()

# Pairwise_Dereplication
A set of tools for dereplication of MAGs via pairwise comparisons

# Requires
* python3 (on flux just module load python-anaconda3)
* networkx 
    * to install on flux:
         ```pip install networkx --user```

# Running the script 
To dereplicate use selectUniqueGenomes.py

This script is responsible for the actual dereplication

Input:
* A tsv file with 4 columns
    * The column names must be bin1,bin2,ANI,Cov
    * The contents of this file is a list of pairwise comparisons and the resulting ANI and Cov values
      frome whatever comparison software was used
* A tsv file from checkm (you can obtain this from using the --tab_table parameter)
* An ANI threshold (default is .99)
* A Coverage threshold (default is .75)
    
Output:
* unique_genomes.txt
  * A list of genomes determined to be unique within the set of comparisons provided
 
 # How it Works
 1. The input table is parsed and all comparisons stored.
 2. The stored comparisons are trimmed. Comparison are only kept if:
     * The comparison is not a genome to itself
     * The ANI value is equal to or above the threshold in BOTH comparisons of the pairwise results
     * The COV value is equal to or above the threshold in BOTH comparisons of the pairwise results
 3. A set is created that has each genome that is present in the remaining comparisons.
 4. From the set, a set of nodes is generated for a graph
 5. An edge is drawn between each node of a comparison that passed trimming
 6. A graph is generated with many nodes not contained in groups of connected nodes, and some groups of connected nodes
 7. If a node has no edges, it is added to a list of unique genomes
 8. If a node has edges, all the genomes connected by edges in the group are compared to each other via a quality metric (completeness -   5*contamination) and the single highest quality genomes is kept from the group and added to the list of unique genomes
 9. The list of Unique genomes is output
   
   
# Generating the pairwise comparions tsv file
If you use pyani_ANIb, the output is a comparison matrix and not a long formatted table.
To get the required file for selectUniqueGenomes.py you can do the following:

Replace with the actual .tab file names

```python convert_table.py ANI.tab COV.tab```


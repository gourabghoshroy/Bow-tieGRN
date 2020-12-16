# Bow-tie Architecture in Gene Regulatory Networks

Author : Gourab Ghosh Roy

Maintainer : Gourab Ghosh Roy <g.ghoshroy@pgr.bham.ac.uk>

<b>Code </b>: 

Requires the networkx 2.x and numpy Python libraries.

There are two sets of Python files. The files in the 'extraction' folder extract GRNs from the source datasets for different species.
For some species there are multiple extraction files. 
The extraction files for GRNs not used in our experiments only print the number of nodes in the extracted network.
The extraction files for GRNs used in our experiments write edges to an output file. 

In the 'bow-tie' folder, the file 'bowtiedecomposition.py' performs the bow-tie decomposition on a GRN. 
The files 'bowtierandomaddedges.py' and 'bowtierandomdeledges.py' respectively add and delete a certain percentage of random edges to the GRN and performs the bow-tie decomposition.
The file 'bowtierandomnetwork.py' performs bow-tie decomposition on random networks of similar node number and degree. 
The output files with number of nodes and regulators in each bow-tie layer are written in the 'output' folder. 

For details and corresponding references see our paper.  


<b>Datasets </b>:
The source GRN datasets are present in the 'datasets' folder.  
The data file for Arabidopsis 'AtRegNet.csv' is not included here because of Github size limitations. Please download it from the AGRIS data server. 
For references of each data source see our paper.


<b>Networks </b>:
The extracted GRNs are included in the 'networks' folder. 
For each GRN edge, the first column is the regulator gene and the second column is the target gene.

 
<b>Output </b>:
The output of the bow-tie decomposition on networks after random edge addition/deletion and on random networks of similar node number and degree are in this folder.
The format is number of nodes in each of the 7 bow-tie layers - CORE, IN, OUT, INTENDRILS, OUTTENDRILS, TUBES, OTHERS, followed by number of regulators in each of the same 7 layers.





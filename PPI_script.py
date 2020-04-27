# Vânia Cecchini
# 21/04/2020
# PPI integration script
# Version for Sébastien

# Import packages
import pandas as pd
import numpy as np
import networkx as nx
import random
import matplotlib.pyplot as plt

# Loading the ppi
ppi = pd.read_csv('hprd.csv',
                    sep='\t',
                    header=0)
print(ppi.head(5))

# Create graph
G = nx.from_pandas_edgelist(ppi, 'Source', 'Target', True, nx.Graph())
print('number of nodes in the network:',G.number_of_nodes())
print('number of edges in the network:',G.number_of_edges())

# Load MD genes
MDgenes = pd.read_csv('MDGenes.csv',
                    sep='\t',
                    header=0)
print(MDgenes.head(5))
print(MDgenes.shape)

# Load miRNA-target interaction file
MTIs = pd.read_csv('miRNAs.csv',
                    sep='\t',
                    header=0)
print(MTIs.head())
print(MTIs.shape)

# Genes have duplicate values as they are targeted by several miRNAs
# We need to remove the duplicate values
# sorting by first name which are the genes (keys)
MTIs.sort_values('genes',
                 inplace = True) # identifying duplicates by the first column
  
# dropping ALL duplicate values
MTIs.drop_duplicates(subset ='genes',
                     keep = False,
                     inplace = True)
print(MTIs.head(5))
print(MTIs.shape)

# Replace all miRNAs names just to yes. At this point, no need to have their names, just the targeting info
# Later we can extract the number of miRNAs that target each gene and add them as weights to the network
MTIs['miRNAs'] = 'Yes'
print(MTIs.head(5))

# Adding node attributes
nx.set_node_attributes(G, MDgenes.set_index('genes').to_dict('index'))
nx.set_node_attributes(G, MTIs.set_index('genes').to_dict('index'))

# Check if the attributes were added
G.nodes(data=True)

# Extracting the different node types according to their attributes
# miRNA Disease nodes
mirna_md_nodes = [n for (n,ty) in \
    nx.get_node_attributes(G,'miRNAs').items() & nx.get_node_attributes(G,'MD').items()]
print('number of mirna disease genes:',len(mirna_md_nodes))

# Disease genes
# First calculate the total number and then subtract the ones that are also mirna disease genes
total_md_nodes = [n for (n,ty) in \
    nx.get_node_attributes(G,'MD').items()]
print('total number of disease proteins:',len(total_md_nodes))
md_nodes = list(set(total_md_nodes) - set(mirna_md_nodes))
print('final number of disease proteins:',len(md_nodes))

# Same for mirna genes
total_mirna_nodes = [n for (n,ty) in \
    nx.get_node_attributes(G,'miRNAs').items()]
print('total number of genes targeted by mirnas:',len(total_mirna_nodes))
mirna_nodes = list(set(total_mirna_nodes) - set(mirna_md_nodes))
print('final number of genes targeted by mirnas:',len(mirna_nodes))

# Normal nodes (remaining nodes)
normal_nodes = list(set(G.nodes()) - set(md_nodes) - set(mirna_nodes) - set(mirna_md_nodes))
print('number of Normal genes in the network:',len(normal_nodes))


# Information on the nodes
print('number of MD genes in the network:',len(md_nodes))
print('number of miRNAs genes in the network:',len(mirna_nodes))
print('number of MD miRNAs genes in the network:',len(mirna_md_nodes))
print('number of Normal genes in the network:',len(normal_nodes))
print('\n')
print('To verify if the calculations are correct:')
total_length = len(md_nodes)+len(mirna_nodes)+len(mirna_md_nodes)+len(normal_nodes)
print('sum of all the obtained nodes:',total_length)
print('number of nodes in the network:',G.number_of_nodes())

# Check if the network is connected
print('is the network connected?:',nx.is_connected(G))

# Draw the network
pos = nx.spring_layout(G)
fig = plt.figure(figsize=(50,50))
nx.draw_networkx(G, pos, nodelist=md_nodes, \
    node_color='red', node_shape='o')
nx.draw_networkx(G, pos, nodelist=mirna_nodes, \
    node_color='yellow', node_shape='o')
nx.draw_networkx(G, pos, nodelist=mirna_md_nodes, \
    node_color='blue', node_shape='o')
nx.draw_networkx(G, pos, nodelist=normal_nodes, \
    node_color='green', node_shape='o')
plt.show()


# Get the degree centrality
# Because the centrality values are not sorted, we can use a function to output the sorted values
import operator

def centrality_sort(centrality_dict):
    return sorted(centrality_dict.items(), key=operator.itemgetter(1), reverse=True)
    

# Which nodes have the highest/lowest degree centrality?
degree_cent = nx.degree_centrality(G) # normalized values
degree_sorted = centrality_sort(degree_cent)
print('5 highest centrality nodes:',degree_sorted[:5]) # 5 highest centrality nodes
print('5 lowest centrality nodes:',degree_sorted[-5:]) # 5 lowest centrality nodes


# Create subnetwork of most connected nodes
# Generate connected components and select the largest:
largest_component = max(nx.connected_components(G), key=len)

# Create a subgraph of G consisting only of this component:
G2 = G.subgraph(largest_component)


print('initial network values:')
print('number of nodes in the network: 9609')
print('number of edges in the network: 39005')
print('\n')
print('connected network values:')
print('number of nodes in the network:',G2.number_of_nodes())
print('number of edges in the network:',G2.number_of_edges())


# Extracting the different node types according to their attributes
# miRNA Disease nodes
mirna_md_nodes = [n for (n,ty) in \
    nx.get_node_attributes(G2,'miRNAs').items() & nx.get_node_attributes(G2,'MD').items()]
print('number of mirna disease genes:',len(mirna_md_nodes))

# Disease genes
# First calculate the total number and then subtract the ones that are also mirna disease genes
total_md_nodes = [n for (n,ty) in \
    nx.get_node_attributes(G2,'MD').items()]
print('total number of disease proteins:',len(total_md_nodes))
md_nodes = list(set(total_md_nodes) - set(mirna_md_nodes))
print('final number of disease proteins:',len(md_nodes))

# Same for mirna genes
total_mirna_nodes = [n for (n,ty) in \
    nx.get_node_attributes(G2,'miRNAs').items()]
print('total number of genes targeted by mirnas:',len(total_mirna_nodes))
mirna_nodes = list(set(total_mirna_nodes) - set(mirna_md_nodes))
print('final number of genes targeted by mirnas:',len(mirna_nodes))

# Normal nodes (remaining nodes)
normal_nodes = list(set(G2.nodes()) - set(md_nodes) - set(mirna_nodes) - set(mirna_md_nodes))
print('number of Normal genes in the network:',len(normal_nodes))


print('number of MD genes in the network:',len(md_nodes))
print('number of miRNAs genes in the network:',len(mirna_nodes))
print('number of MD miRNAs genes in the network:',len(mirna_md_nodes))
print('number of Normal genes in the network:',len(normal_nodes))
print('\n')
print('To verify if the calculations are correct:')
total_length = len(md_nodes)+len(mirna_nodes)+len(mirna_md_nodes)+len(normal_nodes)
print('sum of all the obtained nodes:',total_length)
print('number of nodes in the network:',G2.number_of_nodes())


# Draw the network
pos = nx.spring_layout(G2)
fig = plt.figure(figsize=(50,50))
nx.draw_networkx(G2, pos, nodelist=md_nodes, \
    node_color='red', node_shape='o')
nx.draw_networkx(G2, pos, nodelist=mirna_nodes, \
    node_color='yellow', node_shape='o')
nx.draw_networkx(G2, pos, nodelist=mirna_md_nodes, \
    node_color='blue', node_shape='o')
nx.draw_networkx(G2, pos, nodelist=normal_nodes, \
    node_color='green', node_shape='o')
plt.show()


# Which nodes have the highest/lowest degree centrality?
degree_cent = nx.degree_centrality(G2) # normalized values
degree_sorted = centrality_sort(degree_cent)
print('5 highest centrality nodes:',degree_sorted[:5]) # 5 highest centrality nodes
print('5 lowest centrality nodes:',degree_sorted[-5:]) # 5 lowest centrality nodes


# THIS IS WHERE I AM AT THE MOMENT.
# IN SUM, WE CAN CREATE THE NETWORK, ADD THE NODE ATTRIBUTES AND THEN DIVIDE THE NODES BY TYPES ACCORDING TO THEIR ATTRIBUTES AND PLOT THEM. WE THEN SEE THAT THE NETWORK IS NOT CONNECTED, SO EXTRACT THE MOST CONNECTED COMPONENTS, WE NEED TO GET THE NODE TYPES AGAIN AND THEN WE CAN PLOT OUR SUBNETWORK.
# NOW FOR THE LAST PART, WE NEED TO GET THE DEGREE CENTRALITY FOR EACH TYPE OF NODE. I STARTED WRITING SOME CODE BUT NOTHING THAT WORKS YET. I WON'T PUT IT HERE YET NOT TO CONFUSE YOU BUT LET ME KNOW IF YOU WANT ME TO SEND YOU MY TEST VERSIONS. LET ME KNOW IF YOU HAVE ANY QUESTIONS OR IF SOME PART OF THE CODE IS NOT RUNNING

# TODO : check that the four lists are mutually exclusive

Results = pd.DataFrame(0, index = G2.nodes, columns = ['NN', 'MD', 'RNA', 'MDRNA'])

NeighborhoodSize = 1 #distance

for idx, thisNode in enumerate(G2.nodes):
	L = NeighborhoodSize
	print('node nr %.0f: %s' % (idx, thisNode))
	N = [thisNode]
	while L > 0:
		M = N
		for m in M:
			N = N + list(G2.neighbors(m)) #list of neighbors
		L = L - 1
	N = list(set(N)) #unique values
	for n in N:
		if n in normal_nodes:
			Results.loc[thisNode, 'NN'] = Results.loc[thisNode, 'NN'] + 1
		elif n in md_nodes:
			Results.loc[thisNode, 'MD'] = Results.loc[thisNode, 'MD'] + 1
		elif n in mirna_nodes:
			Results.loc[thisNode, 'RNA'] = Results.loc[thisNode, 'RNA'] + 1
		elif n in mirna_md_nodes:
			Results.loc[thisNode, 'MDRNA'] = Results.loc[thisNode, 'MDRNA'] + 1

print(Results.head())

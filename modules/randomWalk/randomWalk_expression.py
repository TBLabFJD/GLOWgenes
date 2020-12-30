import numpy as np



""" Returns a normalized adjacency matrix from input networks following edge weight transformation proposed by Vanunu et al. 2010

Inputs:

1. Network File:  Path to network file in SIF format - (node A) tab (relationship type) tab (node B) tab (edge weight)
2. Node Filtering: Path to file listing (first column) genes to be considered

"""


def adjacency_matrix(netFile, nodefiltering=None):


	totalInt = []

	# Network to graph

	with open(netFile, "r") as f:
		for line in f:
			linesplit = line.strip().split(" ")
			weight = linesplit[3].strip("}")
			if nodefiltering!=None:
				if(linesplit[0] in nodefiltering and linesplit[1] in nodefiltering):
					totalInt.append("\t".join([linesplit[0],linesplit[1],weight]))
			else:
					totalInt.append("\t".join([linesplit[0],linesplit[1],weight]))



	G = nx.parse_edgelist(totalInt, delimiter="\t", nodetype = str,  data=(('weight',float),))
	

	# Graph to adjacency matrix

	A = nx.adjacency_matrix(G, weight='weight')
	B = A.todense()


	# Normalized weigths
	# B = B / np.amax(B)


	# Normalized adjacency matrix

	total_neighbors = np.sum(B, axis = 0)
	adjacency_matrix = B / np.vectorize(math.sqrt)(np.multiply(total_neighbors, total_neighbors.transpose()))


	return(adjacency_matrix, list(G))








# Adaptation of the Random Walk method implemented by Agrawal et al. 2018

""" Returns list of scores, based on using assoc_gene_vector as the initial vector.
Scores are reported as 0 for the initial seed genes.

Inputs:
adjacency_matrix: adjacency matrix of shape (number_genes, number_genes)
assoc_gene_vector: Numpy array of shape (number_genes,), 1s for disease-associated genes and rest of them representing a tissue association score.
return_prob: The probability of restart in the random walk
"""


def random_walk_scores(adjacency_matrix, assoc_gene_vector, return_prob=0.75):

	ratio = return_prob 
	convergence_metric = 1

	print(float(np.sum(assoc_gene_vector)))
	p0 = assoc_gene_vector/float(np.sum(assoc_gene_vector))
	old_vector = p0
	while (convergence_metric>1e-6):
		new_vector = (1-ratio) * np.dot(adjacency_matrix, old_vector) + ratio * p0
		new_vector = np.squeeze(np.asarray(new_vector))
		convergence_metric = np.linalg.norm(new_vector-old_vector)
		old_vector = np.copy(new_vector)

	assoc_gene_vector_seeds = np.array([ i if i == 1 else 0 for i in assoc_gene_vector ])

	scores = old_vector * (1 - assoc_gene_vector_seeds) 
	
	return scores


#!/usr/bin/env python

########################################################

# GLOWgenes
# Discovery of new disease-gene associations 
# Author: Lorena de la Fuente
# Date: 09/09/2019

########################################################


import random
import argparse
import sys
import csv
import math
import subprocess

import numpy as np
import pandas as pd
import scipy.sparse as sps
from itertools import chain
import networkx as nx
import scipy.stats as stats
import sklearn.metrics as metrics
from distutils.util import strtobool
from collections import defaultdict
from sklearn.preprocessing import minmax_scale

import modules.randomWalk.randomWalk
import modules.prg.prg



## Seed vector construction for propagation

def seedsVector(seeds, nodes, genes_zscore):

	seedvector=[]

	if genes_zscore!=None:
		# adding continues numeric ranking to candidates
		for x in nodes:
			if x in seeds:
				seedvector.append(1)
			elif x in genes_zscore:
				seedvector.append(genes_zscore[x])
			else:
				seedvector.append(0)
	else:
		#binary vector of seeds
		seedvector =  [int(x in seeds) for x in nodes]

	
	seedarray = np.array(seedvector)

	return(seedarray)



## Conversion of gene names to HGNC

def geneSymbolToHGCN(panelgenes):


	## Reading HGNC dicc for conversion

	hgnc_list = []
	hgnc_dicc = {}

	with open("./dictionaries/HGNC_GENES.txt", 'rt', encoding='utf-8') as hgnc:
		for line in hgnc:
			hgnc_list.append(line.strip().split("\t")[1].upper())

	with open("./dictionaries/HGNC_unique_dicc.txt", 'rt', encoding='utf-8') as hgnc:
		for line in hgnc:
			hgnc_dicc[line.strip().split("\t")[1].upper()] = line.strip().split("\t")[0].upper()


	### Panel genes to HGNC

	# genes wo HGNC ids

	panelgenes_hgnc =  []

	for gene in panelgenes:

		gene = gene.upper()
		if gene in hgnc_list:
			panelgenes_hgnc.append(gene)
		elif gene in hgnc_dicc:
			panelgenes_hgnc.append(hgnc_dicc[gene])
		else:
			panelgenes_hgnc.append(gene)


	return(panelgenes_hgnc)




## Computing network parameters

def getNetParameters(graph):


	parameters = [str(graph.number_of_nodes()),
	str(graph.number_of_edges()),
	str(nx.average_clustering(graph))]
	str(nx.number_connected_components(graph)),
	",".join([str(graph.subgraph(c).number_of_nodes()) for c in nx.connected_components(graph)])]


	return(parameters)





## Prediction: network performance metrics

def computePerformance(panelgenes, seeds_list, candidates_list, nodes, adjacencyMatrix, genes_zscore, disGenesOutInt):


	aucDicc = defaultdict(list)
	thresholdDicc = defaultdict(list)
	lst = []	
	dectG = []

	for k in range(0,len(seeds_list)):

		seeds = seeds_list[k]
		testset = candidates_list[k]
		
		## Create numpy array of shape (number_genes,). 
		seedarray = seedsVector(seeds, nodes, genes_zscore)

		# Random Walk with Restart Scoring
		rwwr_scores = modules.randomWalk.randomWalk.random_walk_scores(adjacencyMatrix, seedarray).tolist()
		

		# Prediction vector: 0s for test genes not evaluated (not in network)
		rwwr_scores_incOutInt = [rwwr_scores[x] for x in range(0,len(nodes)) if nodes[x] not in seeds] + [0 for gene in disGenesOutInt if gene not in seeds]
		

		# Labels vector for candidates 
		genes_non_sample = [gene for gene in nodes if gene not in seeds ] + [gene for gene in disGenesOutInt if gene not in seeds ]
		labels_incOut =  [int(gene in testset) for gene in genes_non_sample] 


		# Performance
		lab = np.array(labels_incOut, dtype='int')
		pred = np.array(rwwr_scores_incOutInt, dtype='float')
		fpr, tpr, thresholds = metrics.roc_curve(lab, pred, pos_label=1)
		precision, recall, thresholds = metrics.precision_recall_curve(lab, pred, pos_label=1)

		# auc
		aucDicc["auc"].append(metrics.auc(fpr, tpr))	

		# auprg
		prg_curve = modules.prg.prg.create_prg_curve(lab, pred)
		aucDicc["auprg"].append(modules.prg.prg.calc_auprg(prg_curve))	

		# recall/precision at x-top

		tag_list = ["t", "n", "2",  "4",  "8",  "16",  "32",  "64",  "128",  "256",  "512",  "1024" ]
		ntop_list = [len(testset), len(panelgenes), 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]

		for x in range(len(ntop_list)):

			ntop = ntop_list[x]
			top_n_index = np.argsort(pred)[::-1][:ntop]

			pred_label = [1 if index in top_n_index else 0 for index in range(0,len(pred))]

			f1 = metrics.f1_score(lab, pred_label)
			fbeta = metrics.fbeta_score(lab, pred_label, 2)
			baccuracy = metrics.balanced_accuracy_score(lab, pred_label)
			matthews_corrcoef = metrics.matthews_corrcoef(lab, pred_label)
			recall = metrics.recall_score(lab, pred_label)
			precision = metrics.precision_score(lab, pred_label)

			CM = metrics.confusion_matrix(lab, pred_label)
			tp = CM[1][1]

			lst.append([tag_list[x], ntop, tp, f1, recall, precision, matthews_corrcoef, baccuracy])


			# detected genes

			genesAtN = ",".join([gene for gene in np.array(genes_non_sample)[top_n_index] if gene in testset ])
			dectG.append([tag_list[x], ntop, k, tp, genesAtN])



	cols = ['Type','ntop', '#tp', 'f1', 'recall', 'precision', 'matthews_corrcoef', 'baccuracy']
	dfstats = pd.DataFrame(lst, columns=cols)

	cols = ['Type','ntop', 'iteration', '#tp', 'TPgenes']
	dfdectG = pd.DataFrame(dectG, columns=cols)

	# summary grouped by ntop

	dfstats_mean = dfstats.groupby("Type").mean()
	
	return([np.mean(aucDicc["auc"]), np.mean(aucDicc["auprg"]), dfstats_mean, dfdectG])







## Read disease genes from user input files


def readDiseasegenes(seedfile, temporal, panelapp, trainingratio):

	disease = seedfile.split("/")[-1].split(".txt")[0].split(".tsv")[0]

	diseasegroup = "NA"
	timeprinted = defaultdict()

	if args.timeprinted: # implemented for disgenet   ##### PROBLEMS!!!!!! DIVIDISION WOULD BE OUT OF THE FUNCTION

		panel = pd.read_csv(seedfile,
		sep="\t", header=None).convert_dtypes()

		panel["hgncGene"] = geneSymbolToHGCN(list(panel.iloc[:,0]))

		nseeds = round(panel.shape[0]*trainingratio)

		timeprinted["training"] = panel.iloc[0:nseeds,:].loc[:,'hgncGene'].tolist()
		timeprinted["test"]  = panel.iloc[nseeds:panel.shape[0]+1,:].loc[:,'hgncGene'].tolist()
		timeprinted["time"] = panel.iloc[nseeds,6].tolist()
		panelgenes = panel["hgncGene"]
		

	else:

		if args.panelapp:

			with open(seedfile, "r") as f:
				f.readline()
				panelgenes = [line.strip().split("\t")[0] for line in f if line.strip().split("\t")[1]=="gene"]

			diseasegroup = open(args.input, "r").readlines()[2].strip().split("\t")[6]
	

		else:

			with open(seedfile, "r") as f:
				f.readline()
				panelgenes = [line.strip().split("\t")[0] for line in f]
			
		panelgenes = list(set(geneSymbolToHGCN(panelgenes)))


		return(panelgenes, timeprinted, diseasegroup, disease)




## Main


def main(args):


	# COMMENT CODE

	sys.stdout.write("\nRunning GLOWgenes\n")
	sys.stdout.write("\nInput argument:\n")

	sys.stdout.write(args)

	sys.stdout.write("\nChecking arguments...\n") 



	# reading panel genes 

	sys.stdout.write("\n\n Reading input disease-associated genes ... ")  
	panelgenes, timeprinted, diseasegroup, disease = readDiseasegenes(seedfile = args.input, temporal = args.timeprinted, panelapp = args.panelapp, trainingratio = args.ratio)


	# filtering nodes: for example selection of genes expressed in studied disease tissue
	
	gtofilter = None

	if args.filtering!=None:
		with open(args.filtering, "r") as f:
			gtofilter = [line.strip().split("\t")[1] for line in f]
		gtofilter = geneSymbolToHGCN(gtofilter)


	# ranking candidate genes: for example by tissue expression

	genes_zscore = None

	if args.expnorm!=None:

		gene_symbol = []
		z_score=[]
		z_normalized = []

		with open(args.expnorm, "r") as summary_finalexpression:
			for line in summary_finalexpression:
				gene_symbol.append(line.strip().split("\t")[0])
				z_score.append(float(line.strip().split("\t")[1]))

		z_normalized = minmax_scale(z_score, feature_range=(0, args.cutoff))
		genes_zscore = dict(zip(gene_symbol,z_normalized))


	
	# defining traning/test splits

	testset_list = []
	seeds_list = []
	nseeds = round(len(panelgenes)*args.ratio)


	for r in range(0,20):

		random.shuffle(panelgenes)
		seeds_list.append(panelgenes[0:nseeds])
		testset_list.append(panelgenes[nseeds:len(panelgenes)])



	### Performance + Results for each Network

	netFile = open(args.networks, "r")
	totalstats = pd.DataFrame()
	totalgenes = pd.DataFrame()

	for line in netFile:
		
		path, networkname, category = line.strip().split("\t")

			
		### Create adjacency matrix
		adjacencyMatrix, nodes = modules.randomWalk.randomWalk.adjacency_matrix(netFile=path, nodefiltering=gtofilter)


		### Select disease-genes not in network
		diseasegenesOutInt = list(set(panelgenes) - set(nodes))


		### Check if enough disease genes in network
		seedsByVal = [len(list(set(nodes) & set(seeds_list[k]))) for k in range(0,len(testset_list))]
		


		if all(i >= 1 for i in seedsByVal):

			### Run network performance for cross-validation

			auc, auprg, statsDF, genesDF = computePerformance(panelgenes, seeds_list, testset_list, nodes, adjacencyMatrix, genes_zscore, diseasegenesOutInt)
			statsDF['trainingtype'] = genesDF['trainingtype'] = "Random"
			statsDF['auc'] = auc
			statsDF['auprg'] = auprg


			### Run network performance for time-printed

			if args.timeprinted:

				auc, auprg, timestatsDF, timegenesDF = computePerformance(panelgenes, timeprinted[0], timeprinted[1], nodes, adjacencyMatrix, genes_zscore, diseasegenesOutInt)
				timestatsDF['trainingtype'] = "Temporal"
				timestatsDF['auc'] = auc
				timestatsDF['auprg'] = auprg

				statsDF = pd.concat([statsDF,timestatsDF])


			statsDF['disease'] = genesDF['disease'] = disease
			statsDF['diseasegroup'] = genesDF['diseasegroup'] = diseasegroup
			statsDF['network'] = genesDF['network'] = networkname
			statsDF['knowledgecategory'] = genesDF['knowledgecategory'] = category
			statsDF['#diseasegenes'] = genesDF['#diseasegenes'] = len(panelgenes)
			statsDF['percGenesOutNetwork'] = len(diseasegenesOutInt)/len(panelgenes)*100


			totalstats = totalstats.append(statsDF)
			totalgenes = totalgenes.append(genesDF)

	

			#### Random walk results

			seedarray = seedsVector(panelgenes, nodes, genes_zscore)
			rwwr_results = modules.randomWalk.randomWalk.random_walk_scores(adjacencyMatrix, seedarray).tolist()


			with open(args.output+"/GLOWgenes_"+networkname+"_"+category+".txt", "w") as results:
				ranking=1
				for index in np.argsort(rwwr_results)[::-1]:
					if nodes[index] not in panelgenes:
						results.write("%s\t%f\t%d\n"%(nodes[index],rwwr_results[index],ranking))
						ranking+=1


	performance_threshold = args.output+"/networkEvaluation.txt"
	totalstats.round(decimals=3).to_csv(performance_threshold, sep = "\t", index=True)
	
	detected_genes = args.output+"/detectedgenesEvaluation.txt"
	totalgenes.to_csv(detected_genes, sep = "\t", index=False)


	### Best Network choice 


	bestnets = args.output+"/bestnetsEvaluation.txt"
	totalstats = totalstats.reset_index()
	best = totalstats.loc[totalstats.groupby(["knowledgecategory","trainingtype"])["auprg"].idxmax(), ["network","knowledgecategory","trainingtype"]]
	best.to_csv(bestnets, sep = "\t", index=False)


	### Network specialization and diversity calculation / disease-aware integration

	subprocess.call (["Rscript", "--vanilla", "./modules/diseaseAware_netIntegration.R", "-d", args.output])
	



if __name__ == '__main__':


	#arguments
	parser = argparse.ArgumentParser(description="Prediction of new gene-disease associations with a disease-aware integration of evidence sources")
	parser.add_argument('-i', '--input', help='\t\t Path to disease genes', required=True)
	parser.add_argument('-n', '--networks', help='\t\t Network config file', required=True)
	parser.add_argument('-o', '--output', help='\t\t Output dir', required=True)
	parser.add_argument('-t', '--timeprinted', help='\t\t Training on old genes.', required=False,  action="store_true")
	parser.add_argument('-p', '--panelapp', help='\t\t Disease-associated genes in PanelApp format', required=False,  action="store_true")
	parser.add_argument('-f', '--filtering', help='\t\t Genes to keep', required=False)
	parser.add_argument('-en', '--expnorm', help='\t\t File with expression levels', required=False)
	parser.add_argument('-co', '--cutoff', help='\t\t Cutoff from GTEx cpmlog', required=False, type=float, default='0.01')
	parser.add_argument('-r', '--ratio', help='\t\t Training ratio', required=False, type=float, default='0.7')

	args = parser.parse_args()
	main(args)




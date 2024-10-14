#!/usr/bin/env python

########################################################

# GLOWgenes
# Discovery of new disease-gene associations 
# Author: Lorena de la Fuente

########################################################


import random
import argparse
import sys
import os
import csv
import math
import subprocess
import glob


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
import logging

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
## In this release GLOWgenes accept only human gene names

def geneSymbolToHGNC(panelgenes):


	## Data dir

	pathname_data = os.path.dirname(sys.argv[0])+"/data/"


	## Reading HGNC dicc for conversion

	hgnc_list = []
	hgnc_dicc = {}

	with open(pathname_data+"HGNC_GENES.txt", 'rt', encoding='utf-8') as hgnc:
		for line in hgnc:
			hgnc_list.append(line.strip().split("\t")[1].upper())

	with open(pathname_data+"HGNC_unique_dicc.txt", 'rt', encoding='utf-8') as hgnc:
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
	str(nx.average_clustering(graph)),
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

		tag_list = ["t", "n"] +   [str(x) for x in list(np.geomspace(1, 1024, num=11, dtype=int))]
		ntop_list = [len(testset), len(panelgenes)] + list(np.geomspace(1, 1024, num=11, dtype=int))

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
## GLOWgenes accepts gene lists in column format (plain text) and PanelApp format (also plain text)


def readDiseasegenes(seedfile, temporal, panelapp, trainingratio):

	disease = seedfile.split("/")[-1].split(".txt")[0].split(".tsv")[0]

	diseasegroup = "NA"
	timeprinted = defaultdict()

	if args.timeprinted: # implemented for disgenet

		panel = pd.read_csv(seedfile,
			sep="\t", header=None).convert_dtypes()

		panel["hgncGene"] = geneSymbolToHGNC(list(panel.iloc[:,0]))

		nseeds = round(panel.shape[0]*trainingratio)

		timeprinted["training"] = panel.iloc[0:nseeds,:].loc[:,'hgncGene'].tolist()
		timeprinted["test"]  = panel.iloc[nseeds:panel.shape[0]+1,:].loc[:,'hgncGene'].tolist()
		timeprinted["time"] = panel.iloc[nseeds,6].tolist()
		
		panelgenes = panel["hgncGene"]
		

	else:

		if args.panelapp: # PanelApp format

			with open(seedfile, "r") as f:
				f.readline()
				panelgenes = [line.strip().split("\t")[0] for line in f if line.strip().split("\t")[1]=="gene"]

			diseasegroup = open(args.input, "r").readlines()[2].strip().split("\t")[6]
	

		else:
			with open(seedfile, "r") as f: # Just genes in a column format
				f.readline()
				panelgenes = [line.strip().split("\t")[0] for line in f]
				
		panelgenes = list(set(geneSymbolToHGNC(panelgenes)))

		if(len(panelgenes) == 0):
			sys.stdout.write('ERROR: Seems seeds are not human gene names. Please check documentation\n')
			sys.exit()


	return(panelgenes, timeprinted, diseasegroup, disease)



## Main

def main(args):


	pathname = os.path.dirname(sys.argv[0]) 


	sys.stdout.write("Running GLOWgenes...\n\n")


	sys.stdout.write('Reading input arguments:\n')

	for arg in vars(args):
		arg_value = getattr(args, arg)
		sys.stdout.write("- %s: %s\n" % (arg, arg_value))		

		if arg in ["input","networks", "filtering", "expnorm"] and arg_value!=None:
			if not os.path.isfile(arg_value): 
				 sys.stdout.write("ERROR: '%s' file does not exist\n" %(arg_value))
				 sys.exit()

		if arg == "output":
			if not os.path.isdir(arg_value):
				 sys.stdout.write("ERROR: '%s' output folder does not exist\n" %(arg_value))
				 sys.exit()



	# reading panel genes 

	sys.stdout.write("\n\nReading disease-associated genes...\n")

	try:
		panelgenes, timeprinted, diseasegroup, disease = readDiseasegenes(seedfile = args.input, temporal = args.timeprinted, panelapp = args.panelapp, trainingratio = args.ratio)

		panelgenesFileHGNC=args.output+"/GLOWgenes_trainingGenes.txt"
		with open(panelgenesFileHGNC, "w") as f:
			for ele in panelgenes: 
				f.write(ele+'\n')

	except:
		sys.stdout.write('ERROR: Seeds file was not correct. GLOWgenes accepts either plain text human gene lists in a single column format or PanelApp format file. Please check documentation\n')
		sys.exit()

	sys.stdout.write("- Disease: %s\n- Disease Group: %s\n- #Input disease-associated Genes: %d\n" % (disease, diseasegroup, len(panelgenes)))		
	
	if len(panelgenes) < 40:
		sys.stdout.write('WARNING: Low number of disease-associated genes')
	if len(panelgenes) < 10:
		sys.stdout.write('ERROR: Not enough number of disease-associated genes')
		sys.exit()



	# filtering nodes: for example selection of genes expressed in studied disease tissue
	
	gtofilter = None

	if args.filtering!=None:

		sys.stdout.write("\n\nReading genes to filter networks...\n")

		try:
			with open(args.filtering, "r") as f:
				gtofilter = [line.strip().split("\t")[1] for line in f]
			gtofilter = geneSymbolToHGNC(gtofilter)
			sys.stdout.write("- Genes to keep: %d\n" % len(gtofilter))	
		except:
			sys.stdout.write('ERROR: Gene list without rigth format\n')




	# ranking candidate genes: for example by tissue expression

	genes_zscore = None

	if args.expnorm!=None:

		sys.stdout.write("\n\nReading gene expression values...\n")

		gene_symbol = []
		z_score=[]
		z_normalized = []

		try:
			with open(args.expnorm, "r") as summary_finalexpression:
				for line in summary_finalexpression:
					gene_symbol.append(line.strip().split("\t")[0])
					z_score.append(float(line.strip().split("\t")[1]))

			z_normalized = minmax_scale(z_score, feature_range=(0, args.cutoff))
			genes_zscore = dict(zip(gene_symbol,z_normalized))

		except:
			sys.stdout.write('ERROR: Gene expression file without right format. Please check documentation.\n')

	


	# defining traning/test splits
	
	sys.stdout.write("\n\nRunning GLOWgenes network evaluation and propagation for gene prediction")


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
			statsDF['percgenesinnetwork'] = 100 - (len(diseasegenesOutInt)/len(panelgenes)*100)


			totalstats = totalstats.append(statsDF)
			totalgenes = totalgenes.append(genesDF)

	

			#### Random walk results

			seedarray = seedsVector(panelgenes, nodes, genes_zscore)
			rwwr_results = modules.randomWalk.randomWalk.random_walk_scores(adjacencyMatrix, seedarray).tolist()

			singleNet_output = args.output+"/singleNetworkModeling/"
			try:
				os.mkdir(singleNet_output)
			except OSError:
				pass
			with open(singleNet_output+"GLOWgenes_"+networkname+"_"+category+".txt", "w") as results:
				ranking=1
				for index in np.argsort(rwwr_results)[::-1]:
					if nodes[index] not in panelgenes:
						results.write("%s\t%f\t%d\n"%(nodes[index],rwwr_results[index],ranking))
						ranking+=1


	

	performance_threshold = args.output+"/networkEvaluation.txt"
	totalstats = totalstats.reset_index().round(decimals=3)
	totalstats = totalstats[['disease','diseasegroup','#diseasegenes','network','knowledgecategory', 'percgenesinnetwork',  'trainingtype', 'auprg', 'auc', 'Type', 'ntop','#tp',	'recall',	'precision', 'f1',	'matthews_corrcoef', 'baccuracy']]
	totalstats.to_csv(performance_threshold, sep = "\t", index=False)
	
	detected_genes = args.output+"/detectedgenesEvaluation.txt"
	totalgenes.to_csv(detected_genes, sep = "\t", index=False)


	### Best Network choice 


	bestnets = args.output+"/bestnetsEvaluation.txt"
	best = totalstats.loc[totalstats.groupby(["knowledgecategory","trainingtype"])["auprg"].idxmax(), ["network","knowledgecategory","trainingtype"]]
	best.to_csv(bestnets, sep = "\t", index=False)




	### Network specialization and diversity calculation / disease-aware integration

	sys.stdout.write("\n\nDisease-aware integration of heteregenous networks...\n")
	try:       
		subprocess.call (["Rscript", "--vanilla", pathname+"/modules/diseaseAware_netIntegration.R", "-d", args.output])
	except:
		sys.stdout.write('ERROR: Problems during network integration')
		sys.exit()



	### Variant prioritization for files by in folder (.txt) or specified file

	if args.vcf!=None or args.dir!=None:

		sys.stdout.write("\n\nPrioritizating annotated variants...\n")

		if args.vcf != None:

			if not os.path.isfile(args.vcf): 
				 sys.stdout.write("ERROR: '%s' file does not exist\n" %(arg_value))
				 sys.exit()
			try:       
				subprocess.call (["Rscript", "--vanilla", pathname+"/modules/variantPrioritization.R", "-d", args.output])
			except:
				sys.stdout.write('ERROR: Problems during variant prioritization')
				sys.exit()

		else:

			if not os.path.isdir(args.dir): 
				 sys.stdout.write("ERROR: '%s' directory does not exist\n" %(arg_value))
				 sys.exit()
			else:
				files = glob.glob(args.dir+"/*txt")
				for file in files:
					sys.stdout.write("Analyzing '%s' file...\n" %(file))
					try:       
						subprocess.call (["Rscript", "--vanilla", pathname+"/modules/variantPrioritization.R", "-p", args.output+"/GLOWgenes_prioritization_Random.txt", "-t", panelgenesFileHGNC, "-d", disease, "-v", file, "-a", pathname+"/data"])
					except:
						sys.stdout.write('ERROR: Problems during variant prioritization')
						sys.exit()



	## Intermediary files removal
	os.remove(args.output+"/detectedgenesEvaluation.txt")
	os.remove(args.output+"/networkEvaluation_best.txt")

	sys.stdout.write('GLOWgenes successfully finished!\n')



if __name__ == '__main__':


	#arguments
	parser = argparse.ArgumentParser(description="Prediction of new gene-disease associations by disease-aware integration of evidence sources")
	parser.add_argument('-i', '--input', help='\t\tFile listing already associated disease genes', required=True)
	parser.add_argument('-n', '--networks', help='\t\t Network config file. Three tab-separated fields: network path, network name, network category', required=True)
	parser.add_argument('-o', '--output', help='\t\t Output directory', required=True)
	parser.add_argument('-t', '--timeprinted', help='\t\t Knowledge accumulation approach.', required=False,  action="store_true")
	parser.add_argument('-p', '--panelapp', help='\t\t Disease-associated genes are in PanelApp format', required=False,  action="store_true")
	parser.add_argument('-f', '--filtering', help='\t\t List of candidates list. Edges involving genes not listed here are filtered from networks', required=False)
	parser.add_argument('-en', '--expnorm', help='\t\t Gene expression levels. Two tab-separated fields: gene name, expression level', required=False)
	parser.add_argument('-co', '--cutoff', metavar='float', choices=range(1),  help='\t\t Maximum seed initialization value when considering gene expression levels in the range 0..1', required=False, type=float, default='0.01')
	parser.add_argument('-r', '--ratio', metavar='float', choices=range(1), help='\t\t Training ratio for random training/test splits in the range 0..1', required=False, type=float, default='0.7')
	parser.add_argument('-v', '--vcf', help='\t\t Annotated variants in text format (VEP annotation) to be prioritized by GLOWgenes', required=False)
	parser.add_argument('-d', '--dir', help='\t\t Directory with annotated variants in text format (VEP annotation) to be prioritized by GLOWgenes', required=False)

	args = parser.parse_args()
	main(args)


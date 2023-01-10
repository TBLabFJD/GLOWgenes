# GLOWgenes
Prioritization of gene diseases candidates by disease-aware evaluation of heterogeneous evidence networks
Visit www.glowgenes.org for more information

## Citing
Preprint available at: https://www.medrxiv.org/content/10.1101/2022.10.07.22280759v2

## Requirements

R (tested with version 3.5.0). 
R packages: optparse, caret

Python 2.7 or 3.6

Python packages: numpy (tested with version 1.11.0), pandas (tested with version 0.19.0), scipy (tested with version 0.18.1), sklearn (tested with version 0.0), networkx (tested with version 3.0)

## Obtaining network files
Download network files from: Minguez, Pablo (2022): GLOWgenesNets.zip. figshare. Dataset. https://doi.org/10.6084/m9.figshare.21408393.v1

You could also generate your own networks or selected a subset from theose provided by GLOWgenes

## Editing networks config file
Edit networks_knowledgeCategories.cfg file with your complete directory route to the network files
e.g. substitute PATH by home/pablo/GLOWgenesNets in every line, as in: /PATH/coexpressionCOXPRESdbEXT_HGNCnets.txt

## Running GLOWgenes

usage: GLOWgenes.py [-h] -i INPUT -n NETWORKS -o OUTPUT [-t] [-p]
                    [-f FILTERING] [-en EXPNORM] [-co CUTOFF] [-r RATIO]


python GLOWgenes.py -i diseaseGenes.txt -n networks.cfg -o outputdir -p

Use complete paths to avoid errors

## Parameters

Mandatory parameters:

**-i --input INPUT**
File listing known associated disease genes

**-n --networks NETWORKS**
Evidence network config file. Three tab-separated fields: network path, network name, network category

[Default network config file](networks_knowledgeCategories.cfg)

DEFAULT NETWORK CONFIG FILE IS LOCATED AT TEST FOLDER

**-o --output OUTPUT**
Output directory

**-p, --panelapp**       
Disease-associated genes in PanelApp format

Gene Panels from PanelApp can be download from https://panelapp.genomicsengland.co.uk/panels/.

**-t, --timeprinted**     
Knowledge accumulation approach.
  
**-f FILTERING, --filtering FILTERING**
List of candidate genes. Edges involving genes not listed here are filtered from networks
  
**-en EXPNORM, --expnorm EXPNORM**
Expression levels file. Two tab-separated fields: gene name, expression level               


**-co CUTOFF, --cutoff CUTOFF**
Maximum seed initialization value when considering gene expression levels. Range 0-1
  
  
**-r RATIO, --ratio RATIO**
Training ratio for random training/test splits

## Running an example
Within directory example you have full intructions to test GLOWgenes

############################################################################################
## Disease-aware integration of funcional knowledge for the discovery of new disease-genes
## Author: Lorena de la Fuente
## 28-12-2020
############################################################################################

library(optparse)
library(caret)


# Network Specialization & Diversity 

calculateSpecScore <- function(detectedgenes){
  
  agene=table(unlist(detectedgenes))
  mymatrix = matrix(nrow = length(detectedgenes), ncol = length(agene))
  colnames(mymatrix) = names(agene)
  rownames(mymatrix) = names(detectedgenes)
  
  for (x in names(detectedgenes)){
    mymatrix[x,] <- as.numeric(colnames(mymatrix) %in% detectedgenes[[x]])
  }
  
  mymatrix_DIV = mymatrix/rowSums(mymatrix)
  
  mymatrix_DIV[is.na(mymatrix_DIV)] <- 0
  
  Hj = apply(mymatrix_DIV,1,function(x) -1*sum(x*log2(x), na.rm=TRUE)) # Diversity of the transcriptome of each tissue
  
  pi = apply(mymatrix_DIV, 2, mean)
  
  Si = 1/nrow(mymatrix_DIV) * colSums(t(t(mymatrix_DIV)/pi)*log2(t(t(mymatrix_DIV)/pi)),na.rm=TRUE)
  
  gj = rowSums(t(t(mymatrix_DIV)*Si)) / log2(length(detectedgenes))

  return(list(Hj, gj))
  
}


# Network Specialization & Diversity 

kcIntegration <- function(metrics_best, dir){
  
  
  # Reading RWWR results
  rwwr.list = sapply(unique(metrics_best$network), function(y) list.files(path=dir, pattern = paste("GLOWgenes_",y, sep = ""), full.names = T))
  rwwr.results = NULL
  genes = NULL
  for (n in 1:length(rwwr.list)){
    
    net = names(rwwr.list)[n]
    rwwr.results[[net]] = read.delim(file=rwwr.list[n], header = FALSE, row.names = 1)
    genes = unique(c(genes,rownames(rwwr.results[[net]] )))
    
  }
  
  
  # RWWR matrix
  myresults=NULL
  myranking = NULL
  
  for ( element in rwwr.results ){
    
    myresults=cbind(myresults, element[genes,"V2"])
    myranking=cbind(myranking, element[genes,"V3"])
    
  }
  
  rownames(myranking) = rownames(myresults) = genes
  colnames(myranking) = colnames(myresults) = names(rwwr.results)
  
  
  # KCs integration
  
  metrics_best.n = metrics_best[metrics_best$Type=="n",]
  rownames(metrics_best.n) = metrics_best.n$network
  
  scoreNorm = scale(myresults, center = T, scale = T)
  scoreNorm_imputed <- apply(scoreNorm, 2, function(x)  {x [ is.na(x)] <-  min(x, na.rm=T); return(x)})
  weightedScoreNorm = t(t(scoreNorm_imputed) * as.numeric(metrics_best.n[colnames(scoreNorm_imputed),"specificity"]) * as.numeric(metrics_best.n[colnames(scoreNorm_imputed),"recall"]))
  weigthedFinalNormImp = data.frame(gene = rownames(weightedScoreNorm), combinedWeigthed=apply(weightedScoreNorm, 1, mean, na.rm = TRUE))
  weigthedFinalNormImp = weigthedFinalNormImp[order(weigthedFinalNormImp$combinedWeigthed, decreasing = TRUE),]
  weigthedFinalNormImp$rank = 1:nrow(weigthedFinalNormImp)

  return(weigthedFinalNormImp)  
  
}



##### Main


option_list=list(
  make_option(c("-d","--dir"), type="character", help="Results directory")
)

opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args



# Loading data

mydetgenes = read.delim2( file = paste(opt$dir,"detectedgenesEvaluation.txt",sep="/"), stringsAsFactors = F)
metrics = read.delim2(file = paste(opt$dir,"networkEvaluation.txt",sep="/"), stringsAsFactors = F)





# KC integration by strategy

best = read.delim2(file = paste(opt$dir,"bestnetsEvaluation.txt",sep="/"), stringsAsFactors = F)
best.bystrategy = split(best$network, best$trainingtype)


metrics.best.spec.bystrategy = lapply(names(best.bystrategy), function(x) {
  
  mydetgenes.best = mydetgenes[ mydetgenes$trainingtype==x & mydetgenes$network%in%best.bystrategy[[x]],]
  metrics.best = metrics[ metrics$trainingtype==x & metrics$network%in%best.bystrategy[[x]],]
  
  specificity.best = NULL
  
  for (ntop in unique(mydetgenes.best$Type)){
    
    genes.topn = mydetgenes.best[mydetgenes.best$Type==ntop,]
    genes.topn.list = split(x = genes.topn, f = genes.topn$iteration)
    
    # for each random_sampling
    
    spec = lapply(genes.topn.list, function(z) { 
      
      genelist = sapply(z$TPgenes, function(y) strsplit(y, split = ",", fixed = T))
      
      names(genelist) = z$knowledgecategory
      
      spec.ite = calculateSpecScore(detectedgenes = genelist)
      
      return(spec.ite)
      
    }
    )
    
    gj_mean = rowMeans(sapply(spec, function(m) m[[2]]))
    hj_mean = rowMeans(sapply(spec, function(m) m[[1]]))
    
    specificity.best = rbind(specificity.best, data.frame(knowledgecategory = names(gj_mean), diversity = hj_mean, specificity = gj_mean, Type=ntop, trainingtype=x))
    
  }
  
  metrics.best.spec = merge(metrics.best, specificity.best, by=c("Type","knowledgecategory")) 
  
  
  # Load random-walk results & Disease-aware Integration 
  
  integration.results = kcIntegration(metrics_best = metrics.best.spec, dir = opt$dir)
  write.table(integration.results, file=paste(opt$dir,"/GLOWgenes_prioritization_", x, ".txt", sep = "") , row.names = FALSE, col.names = FALSE, quote = FALSE,sep = "\t")
  
  
  return(metrics.best.spec)
  
  })






# Write Best Network Performance

metrics.best.spec =  do.call("rbind", metrics.best.spec.bystrategy)
write.table(metrics.best.spec, file=paste(opt$dir,"networkEvaluation_best.txt", sep = "/"),  row.names = F, quote = FALSE,sep = "\t")









##############################################################
# Variant prioritization by GLOWgenes candidate prediction 
# Author: Lorena de la Fuente
# Date: 09-05-20
###############################################################

library(optparse)


##### Functions


# HGNC - all genes to HGNC  

tohgnc <- function(list, dir){

  # Reading information
  hgnc = read.delim2(paste(dir,"/","HGNC_unique_dicc.txt", sep=""), stringsAsFactors = F, header = F)
  hgncgenes = read.delim2(paste(dir,"/","HGNC_GENES.txt", sep=""), stringsAsFactors = F, header = T)
  hgncgenes$upper = toupper(hgncgenes$Approved.symbol)
  hgnc = hgnc[!is.na(hgnc$V2),]
  rownames(hgnc) <- hgnc$V2
  hgnc$upper <- toupper(hgnc$V2)
  
  approved = hgncgenes[match(x = toupper(list), table = hgncgenes$upper), "Approved.symbol"]
  approved[is.na(approved)] = hgnc[match(x = toupper(list)[is.na(approved)], table = hgnc$upper), "V1"]
  return(approved)
  
}



##### Main

option_list=list(
  make_option(c("-p","--prioritization"), type="character", help="Prioritization file"),
  make_option(c("-t","--training"), type="character", help="Training genes"),
  make_option(c("-d","--disease"), type="character", help="Disease name/id"),
  make_option(c("-v","--variants"), type="character", help="Annotated variants in text format"),
  make_option(c("-a","--data"), type="character", help="Data directory where HGNC dictionary"),
  make_option(c("-o","--output"), type="character", help="Output directory. Must be created"),
  make_option(c("-f","--filter"), type="character", help="Filtering of variants by pathogenic", default=TRUE))

opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args


if (is.null(opt$output)){
  opt$output=dirname(opt$variants)
}else{if (!dir.exists(opt$output)){cat("ERROR: Output directory does not exist")}}


## Loading GLOWgenes Prioritization Results
glowgenes.ranking = read.delim(opt$prioritization, header=F, stringsAsFactors = F)

## Loading training genes
training.genes = read.delim(opt$training, header=F, stringsAsFactors = F)[,1]

## Loading variants
myvariants = read.delim(file=opt$variants, header = T, sep = "\t", stringsAsFactors = F)
basename = tools::file_path_sans_ext(basename(opt$variants))

# Prioritizating variants by GLOWgenes ranking
myvariants$SYMBOLOK = tohgnc(list = myvariants$SYMBOL, dir = opt$data)
nameprior = "GLOWgenesRanking"
myvariants[,nameprior] <- match(myvariants$SYMBOLOK, glowgenes.ranking$V1)  
myvariants[myvariants$SYMBOLOK%in%training.genes, nameprior] <- 0

# Writing results
output=paste(opt$output, "/", basename, "_GGPrioritized", ".txt", sep="")
write.table(myvariants, file=output, sep = "\t", quote=F, row.names = F) 



# Prioritization plus variant filtering

if(opt$filter==TRUE){
  
  #* MAX_AF <0.05 o vacias OR gnomad popmax  <0.05 o vacias
  #* CADD PHRED >10 o vacias
  #* CANONICAL == YES
  #* GLOWgenes prioritization in any panel < 500
  
  myvariants.filtered <- subset(myvariants, 
                             (gnomADg_AF_POPMAX <= 0.05 | is.na(gnomADg_AF_POPMAX) |  MAX_AF <= 0.05 | is.na(MAX_AF)) &
                               (CADD_PHRED >= 10 | is.na(CADD_PHRED)) &
                               CANONICAL == "YES" & GLOWgenesRanking<500)
  myvariants.filtered <- myvariants.filtered[!is.na(myvariants.filtered$SYMBOLOK),]
  myvariants.filtered$SYMBOLOK <- NULL
  

  outputfiltered=paste(opt$output, "/", basename, "_GGPrioritizedFiltered", ".txt", sep="")
  write.table(myvariants.filtered, file=outputfiltered, sep = "\t", quote=F, row.names = F) 

}


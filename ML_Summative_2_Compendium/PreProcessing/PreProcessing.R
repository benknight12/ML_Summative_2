## Load packages and data as required
### This script aims to do initial scaling and preprocessing of data
## Load packages
library(data.table)
library(dplyr)

read.dataset=function(filename, ...) {
  cat("reading", basename(filename), "... ")
  ## read in tab-delimited spreadsheet
  x=fread(
    filename,
    header=T,
    stringsAsFactors=F,
    sep="\t",
    check.names=F,
    ...)
  ## remove any duplicate rows (identified by the first column)
  x=x[match(unique(x[[1]]), x[[1]]),]
  ## make the first column the rownames of the data frame
  x=data.frame(x,row.names=1,stringsAsFactors=F,check.names=F)
  cat(nrow(x), "x", ncol(x), "\n")
  x
}


minclinical.dat<-as.data.frame(read.dataset("../rawdata/clinical.txt"))
table(is.na(minclinical.dat))
## 884 missing values
minmethylation.dat<-t(as.data.frame(read.dataset("../rawdata/methylation.txt")))
table(is.na(minmethylation.dat))
## 10696 missing values
minmirna.dat<-t(as.data.frame(read.dataset("../rawdata/mirna.txt")))
table(is.na(minmirna.dat))
## 0 missing values
minmrna.dat<-t(as.data.frame(read.dataset("../rawdata/mrna.txt")))
table(is.na(minmrna.dat))
## 0 missing values
minmutations.dat<-t(as.data.frame(read.dataset("../rawdata/mutations.txt")))
table(is.na(minmutations.dat))
## 0 missing values
minprotein.dat<-t(as.data.frame(read.dataset("../rawdata/protein.txt")))
table(is.na(minprotein.dat))
## 884 missing values
## Not sure what cnv acc is
mincnv.dat<-(as.data.frame(read.dataset("../rawdata/cnv.txt")))
table(is.na(mincnv.dat))
## no missing values
## File which links genes to proteins
linkdat<-as.data.frame(read.dataset("../rawdata/protein-annotation.txt"))
print("6")
dim(minmethylation.dat)
dim(minmirna.dat)
dim(minmrna.dat)
dim(minmutations.dat)
dim(minprotein.dat)
clinical <- minclinical.dat

omics <- list(as.data.table(minmethylation.dat), 
            as.data.table(minmirna.dat),
            as.data.table(minmrna.dat), 
            as.data.table(minmutations.dat),
            as.data.table(minprotein.dat))

names(omics) <- c("methylation",
                "mirna","mrna","mutations",
                "protein")

## replace rownames again
rownames(omics$methylation)<-rownames(minmethylation.dat)
rownames(omics$mirna)<-rownames(minmirna.dat)
rownames(omics$mrna)<-rownames(minmrna.dat)
rownames(omics$mutations)<-rownames(minmutations.dat)
rownames(omics$protein)<-rownames(minprotein.dat)
print("7")
## Scaling functions
minmax <- function(x,...){
  scaled<- (x-min(x,...))/(max(x,...)-min(x,...))
}
scale_na <-function(x,...){
  if(!is.logical(x)) minmax(x,...)
  else return(x)
}
### Do some data transforms on all omics
for (i in 1:length(omics)){
  ### Assume all features numeric unless binary 
  is.binary = sapply(omics[[i]], function(v) all(v==0|v==1,na.rm=T))
  binary.cols = colnames(omics[[i]])[is.binary]
  omics[[i]][,(binary.cols) := lapply(.SD,as.logical), .SDcols = binary.cols]
  ### Remove data with more than 10% missing data
  p.missing = sapply(omics[[i]], function(v) mean(is.na(v)))
  missing.cols = colnames(omics[[i]])[p.missing > 0.10]
  omics[[i]][,(missing.cols) := NULL]
  ## Scale the omic sets
  lapply(omics[[i]], function(x) scale_na(x, na.rm=T))
  print(dim(omics[[i]]))
}
## Produce the following datasets:
pfivalues<- as.data.frame(cbind(pfi=clinical$pfi, names=rownames(clinical)))
pfitime<- as.data.frame(cbind(pfi.time=clinical$pfi.time, names=rownames(clinical)))
## Extract participant id from longer format
get.id <- function(id){
  sub("TCGA-[^-]+-([^-]+)-.*", "\\1", id)
}



## Add the pfis to the data
addpfi <- function(x){
  # Find ids present in both clinical and geneomic data so we can merge them
  similar = intersect(pfivalues$names, get.id(rownames(x)))
  # Create the shrunk subsets to just the shared ids
  x=x[match(similar, get.id(rownames(x))),]
  pfis=pfivalues[match(similar, pfivalues$names),]
  # Merge them
  link.dat=data.frame(cbind(x, pfis))
  rownames(link.dat) <- link.dat$names
  return(as.data.frame(link.dat[,1:ncol(link.dat)-1]))
}
print("8")
omics2<-lapply(omics, addpfi)
## This maybe reduced size a bit not sure
minmethylation.dat<-omics2$methylation
minmirna.dat<-omics2$mirna
minmrna.dat<-omics2$mrna
minmutations.dat<-omics2$mutations
minprotein.dat<-omics2$protein


#### Clinical preprocessing

## Section to dummy include categorical clinical variables
fact_columns<-names(clinical[,lapply(clinical,class)=="character"])
# must include nas as specific category which we will impute later
clinical[,fact_columns] <- lapply(clinical[,fact_columns], function(col) {factor(addNA(col), exclude = NULL)})
# model.matrix used for dummy variable creation
dummy_df <- as.data.frame(model.matrix(~ . - 1, data = clinical[, fact_columns]))
#ensure unique and no space column names
colnames(dummy_df) <- make.names(gsub(" ", "_", colnames(dummy_df)), unique = TRUE)
# add back to clinical
dummy_df<-lapply(dummy_df,as.logical)
clinical <- cbind(clinical[, lapply(clinical,class)!="factor"], dummy_df)


saveRDS(minmethylation.dat,"../preprocessed/methylation.Rda")
saveRDS(minmrna.dat,"../preprocessed/minmrna.Rda")
saveRDS(minmirna.dat,"../preprocessed/minmirna.Rda")
saveRDS(minmutations.dat,"../preprocessed/mutations.Rda")
saveRDS(minprotein.dat,"../preprocessed/protein.Rda")
saveRDS(clinical,"../preprocessed/clinical.Rda")
saveRDS(linkdat,"../preprocessed/linkdat.Rda")

## Add the pfitimess to the data
addpfitime <- function(x){
  # Find ids present in both clinical and geneomic data so we can merge them
  similar = intersect(pfivalues$names, get.id(rownames(x)))
  # Create the shrunk subsets to just the shared ids
  x=x[match(similar, get.id(rownames(x))),]
  pfitimes=pfitime[match(similar, pfitime$names),]
  # Merge them
  link.dat=data.frame(cbind(x, pfitimes))
  rownames(link.dat)<-link.dat$names
  return(as.data.frame(link.dat[,1:ncol(link.dat)-1]))
}
print("9")
omics3<-lapply(omics, addpfitime)

tminmethylation.dat<-omics3$methylation
tminmirna.dat<-omics3$mirna
tminmrna.dat<-omics3$mrna
tminmutations.dat<-omics3$mutations
tminprotein.dat<-omics3$protein

saveRDS(tminmethylation.dat,"../preprocessed/methylation_time.Rda")
saveRDS(tminmrna.dat,"../preprocessed/minmrna_time.Rda")
saveRDS(tminmirna.dat,"../preprocessed/minmirna_time.Rda")
saveRDS(tminmutations.dat,"../preprocessed/mutations_time.Rda")
saveRDS(tminprotein.dat,"../preprocessed/protein_time.Rda")

## Linking mutation (genes acc) to proteins
rownames(minmethylation.dat)
print("10")
## Derive the common linked 'features'
minmrna.dat$names <- rownames(minmrna.dat)
minprotein.dat$names <- rownames(minprotein.dat)
genanncommon.features=intersect(names(minmrna.dat), c(linkdat$gene,"names"))
genprotcommon.features=intersect(genanncommon.features, names(minprotein.dat))
genreduced <- minmrna.dat[,genprotcommon.features]
proreduced <-minprotein.dat[,genprotcommon.features]


pfid<-addpfi(genreduced)
compiled <- merge(pfid, proreduced,by="names")
rownames(compiled)<- compiled$names
compiled <- compiled[,-1]


## Save the linked dataset
saveRDS(compiled,"../preprocessed/gene_prot_joined_set.Rda")




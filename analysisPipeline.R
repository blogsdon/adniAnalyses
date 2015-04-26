#get coordinates
source('reduceToFunctionalAnnotation2.R')
for (i in 22:1){
	reduceToFunctionalAnnotation2(i)
}



#Get gene lists
require(synapseClient)
synapseLogin()


geneListSyn <- synGet('syn2673273')

geneList <- read.table(geneListSyn@filePath,sep=',',header=T)

geneListFinal <- unique(toupper(as.character(geneList$Gene)))

geneListFinal2 <- unique(sapply(geneListFinal,function(x) strsplit(x,split=' ')[[1]][1]))

geneListFinal3 <- geneListFinal2[!is.na(geneListFinal2)]

corfs<-intersect(grep(pattern='ORF',geneListFinal3),grep(pattern='C',geneListFinal3))

fixCorfs <- function(x){
	a <- strsplit(x,split='ORF')[[1]]
	b <- paste(a[1],'orf',a[2],sep='')
	return(b)
}

isValidAlias <- function(gene){
	require(org.Hs.eg.db)
	#convert alias to entrez gene id
	xx <- as.list(org.Hs.egALIAS2EG)
	xx <- xx[!is.na(xx)]
	egid <- xx[[gene]][1]



	y <- org.Hs.egCHR
	# Get the entrez gene identifiers that are mapped to a chromosome
	mapped_genes <- mappedkeys(y)
	# Convert to a list
	yy <- as.list(y[mapped_genes])

	chr <- yy[[egid]][1]

	return(list(egid=egid,chr=chr))

}


geneListFinal3[corfs] <- sapply(geneListFinal3[corfs],fixCorfs)

geneListFinal3[geneListFinal3=='1-SEP'] <- 'SEPT1'

geneTable <- matrix(NA,length(geneListFinal3),3)
colnames(geneTable) <- c('gene','eg','chr')
geneTable[,1] <- geneListFinal3;

for(i in 1:length(geneListFinal3)){
	res <- NA
	print(i)
	try(res <- isValidAlias(geneListFinal3[i]),silent=T)
	if(!is.na(res)){
		geneTable[i,'eg'] <- res$egid
		geneTable[i,'chr'] <- res$chr
	}
}

save(geneTable,file='genetable.rda')
geneTable2 <- geneTable[order(as.numeric(geneTable[,3])),]

#populate gene table

#filter gene lists for eg ids
rareVariantData <- vector('list',nrow(geneTable2))
names(rareVariantData) <- geneTable2[,1]

count <- 1;
currentChr <- as.numeric(geneTable2[count,3])
for (i in 1:22){
	snpFunction <- read.table(paste('adniWGSchr',i,'functional.vout',sep=''))
	while(currentChr==i){
	currentChr <- as.numeric(geneTable2[count,3])		
		try(rareVariantData[[count]] <- extractSnps(geneTable2[count,1],snpFunction),silent=TRUE)
		count <- count+1
		print(count)
	}
}

save(rareVariantData,file='rareVariants.rda')



#run plink to extract data 
a <- extractSnps()




#filter snps down to only functional snps
#b <- functSnps()

#run analysis
res <- runAnalysis()



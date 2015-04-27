#extractSnps.R
extractSnps2 <- function(gene,snpFunction){
	require(org.Hs.eg.db)
	#convert alias to entrez gene id
	xx <- as.list(org.Hs.egALIAS2EG)
	xx <- xx[!is.na(xx)]
	egid <- xx[[gene]][1]
	if(length(egid)>0){
		#test if there are any function snps
		#in gene
		n.snp <- sum(snpFunction$Gene==egid)
		if(n.snp>0){
			cat('There are',n.snp,'functional snps in',gene,'\n')
		}else{
			stop('No functional snps in gene',gene,'\n')
		}
		snpLocation <- as.character(snpFunction[snpFunction$Gene==egid,'Location'])
		snpLocation <- strsplit(snpLocation,'\\:')

		chr <- snpLocation[[1]][1];
		positions <- sapply(snpLocation,function(x) return(x[2]))
    position1 <- min(positions);
    position2 <- max(positions);
    cat(chr,position1,position2,'R1\n',file=paste0(gene,'.range'))
		#range <- cbind(rep(chr,n.snp),positions,positions,paste('R',1:n.snp,sep=''))
		#write(t(range),file=paste(gene,'.range',sep=''),ncolumns=4)

		#extract snps from relevant plink file
		str <- paste('plink --bfile adniWGSchr',chr,' --extract ',paste(gene,'.range',sep=''),' --range --recodeA --out ',gene,' --noweb',sep='')
		system(str)
		res <- read.table(paste(gene,'.raw',sep=''),header=T)
		rownames(res) <- res[,1];
		res <- res[,-c(1:6)];


		str2 <- paste('rm ',gene,'.log',sep='')
		cat(str2,'\n')
		system(str2)

		str2 <- paste('rm ',gene,'.nosex',sep='')
		cat(str2,'\n')
		system(str2)

		str2 <- paste('rm ',gene,'.range',sep='')
		cat(str2,'\n')
		system(str2)

		str2 <- paste('rm ',gene,'.raw',sep='')
		cat(str2,'\n')
		system(str2)


		return(res)
	}else{
		stop('not a valid gene symbol')
	}
}




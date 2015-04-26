#reduceToFunctionalAnnotation.R
reduceToFunctionalAnnotation2 <- function(chr){
	x1 <- read.delim(paste('adniWGSchr',chr,'.vout',sep=''),header=T,skip = 8)
	conse <- x1$Consequence
	func <- union(grep('splice_donor',conse),grep('stop',conse))
	func <- union(func,grep('splice_acceptor',conse))
	func <- union(func,grep('missense',conse))
	func <- union(func,grep('intron',conse))
	func <- union(func,grep('synonymous',conse))
	func <- union(func,grep('splice',conse))
	x2 <- x1[func,]
	x3 <- x2[!duplicated(x2$Location),]
	write.table(x3,file=paste('adniWGSchr',chr,'functional2.vout',sep=''))
}


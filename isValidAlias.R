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
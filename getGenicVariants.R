###get all genic variants
require(dplyr)

for (i in 1:22){
  #load .vouts
  x <- read.table(paste0('adniWGSchr',i,'functional2.vout',row.names=NULL),stringsAsFactors=F)
    
  #exclude intron_variant
  y <- dplyr::filter(x,Consequence!='intron_variant'&Consequence!='intron_variant,nc_transcript_variant')
  #save positions to file
  internal <- function(x,b){
    z <- strsplit(x,':')[[1]]
    z <- c(z[1],z[2],z[2],b)
    return(z)
  }
  for(j in 1:nrow(y)){
    cat(paste(internal(y$Location[j],paste0('chr',i,'_',j)),collapse=' '),'\n',file='testTest.out',append=TRUE,sep='')
  }
}

str <- 'plink --bfile adniWGS --extract testTest.out --range --recodeA --out all_var --noweb'
system(str)


 
#' @title Filterin
#' @description This function filter matrix raw count 
#' @param minDelta filtering value for the minimun difference value between genes WMT&rib and genes without MT and without RB
#' @param minNGene filtering value for the minimun number of gene
#' @param original matrix name without annotation 
#' @param menoRibo matrix name with annotation
#' @param separator matrix separator value 
#' @param format matrix format 
#' @author Luca Alessandri, alessandri.luca [at] gmail [dot] com, University of Torino
#' @return filtered matrix table
#' @examples
#' \dontrun{
#'    
#' }
#'
#' @export

deltaFilter=function(minDelta,minNGene,original,menoRibo,separator,format){

mainMatrix=read.table(paste(original,".",format,sep=""),sep=separator,header=TRUE,row.names=1)
menoRibo=read.table(paste(menoRibo,".",format,sep=""),sep=separator,header=TRUE,row.names=1)

    
xCoord=log10(colSums(menoRibo))
b=menoRibo
b[b<3]=0
b[b>=3]=1
yCoord=colSums(b)


xCoord2=log10(colSums(mainMatrix))
b2=mainMatrix
b2[b2<3]=0
b2[b2>=3]=1
yCoord2=colSums(b2)


cells=intersect(which(yCoord2>=minNGene),which((yCoord2-yCoord)>minDelta))
pdf("geneIntersection.pdf")

 plot(yCoord2,yCoord2-yCoord,cex=0.2,pch=19,col="purple",xlab="# of genes", ylab="genesWMT&rib - genes-MT-RB")
 points(yCoord2[cells],(yCoord2-yCoord)[cells],cex=0.2,pch=19,col="red")


  
 dev.off()

 #////////////////
 
 
 write.csv(menoRibo[,cells],"filtered_annotated_lorenz_naive_penta2_0.csv")

}

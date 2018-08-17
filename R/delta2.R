
#' @title Fi
#' @description This function filter matrix raw count
#' @param threshold value to estimate if a gene is significatively expressed
#' @param minDelta filtering value for the minimun difference value between genes WMT&rib and genes without MT and without RB
#' @param minNGene filtering value for the minimun number of gene
#' @param original matrix name without annotation
#' @param menoRibo matrix name with annotation
#' @param separator matrix separator value
#' @param format matrix format
#' @param wf ,if this parameter is setted to 0 your filtering will be minDelta and minNGene based, otherwise it will take the top wf cells with the higher number of genes significatively expressed 
#' @author Luca Alessandri, alessandri.luca [at] gmail [dot] com, University of Torino
#' @return filtered matrix table
#' @examples
#' \dontrun{
#'    TOO BE MADE
#' }
#'
#' @export

deltaFilter2 <- function(threshold,minDelta,minNGene,original,menoRibo,separator,format,wf){
menoRiboName=menoRibo
mainMatrix=read.table(paste(original,".",format,sep=""),sep=separator,header=TRUE,row.names=1)
menoRibo=read.table(paste(menoRibo,".",format,sep=""),sep=separator,header=TRUE,row.names=1)

if(wf==0){    
b=menoRibo
b[b<threshold]=0
b[b>=threshold]=1
yCoord=colSums(b)


#xCoord2=log10(colSums(mainMatrix))
b2=mainMatrix
b2[b2<threshold]=0
b2[b2>=threshold]=1
yCoord2=colSums(b2)


cells=intersect(which(yCoord>=minNGene),which((yCoord2-yCoord)>minDelta))
pdf(paste(menoRiboName,"_geneIntersection.pdf",sep=""))

 plot(yCoord,yCoord2-yCoord,cex=0.2,pch=19,col="purple",xlab="# of genes", ylab="genesWMT&rib - genes-MT-RB",main=paste(ncol(menoRibo[,cells])," Cells taken"))
 points(yCoord[cells],(yCoord2-yCoord)[cells],cex=0.2,pch=19,col="red")


  
 dev.off()

 #////////////////

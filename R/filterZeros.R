#' @title Plotting the distribution of zeros in cells eliminating all genes without counts
#' @description This function plots the zeros distributions in cells and removes genes without counts
#' @param file, a character string indicating the path of the file tab delimited  of cells un-normalized expression counts
#' @param threshold, a number from 0 to 1 indicating the fraction of max accepted zeros in each gene. 0 is set as default and it eliminates only genes having no expression in any cell.
#' @param sep, separator used in count file, e.g. '\\t', ','
#' @return a PDF providing zeros distributions before removal of all genes without counts. A tab delimited file with the prefix *filtered* in which the filtered data are saved.
#' @examples
#' \dontrun{
#'     #downloading fastq files
#'     system("wget http://130.192.119.59/public/singlecells_counts.txt.gz")
#'     system("gzip -d singlecells_counts.txt.gz")
#'     filterZeros(file=paste(getwd(),"singlecells_counts.txt",sep="/"), threshold=0.1, sep="\t")
#' }
#' @export
filterZeros <- function(file, threshold=0, sep){
data.folder=dirname(file)
  setwd(data.folder)
positions=length(strsplit(basename(file),"\\.")[[1]])
matrixNameC=strsplit(basename(file),"\\.")[[1]]
counts.table=paste(matrixNameC[seq(1,positions-1)],collapse="")
  matrixName=counts.table
file.type=strsplit(basename(basename(file)),"\\.")[[1]][positions]

counts.matrix=matrixName
  counts <- read.table(paste(counts.matrix,".",file.type,sep=""), sep=sep, header=T, row.names = 1)
  counts.sum <- apply(counts, 1, function(x){
    length(which(x > 0))
  })
  max.zeros <- dim(counts)[2]*threshold
  counts.n0 <- counts[which(counts.sum > max.zeros),]
  counts.sum0 <- apply(counts.n0, 1, function(x){
    length(which(x > 0))
  })

  cat("\n",paste("Out of ", dim(counts)[1]," genes ",dim(counts.n0)[1]," are left after removing genes with no counts", sep=""),"\n")
  pdf(paste("Non-zeros_distribution_",sub(".txt$","",counts.matrix),".pdf",sep=""))
      hist(counts.sum, col=rgb(1,0,0,0.5), main="",  xlab="# of cells", ylab="Frequency of genes with more than 0 counts")
      hist(counts.sum0, col=rgb(0,0,1,0.5), add=T)
      box()
  dev.off()
  write.table(counts.n0, paste("filtered_",counts.matrix,".",file.type,sep=""), sep="\t", col.names=NA)

}

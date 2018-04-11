#' @title Plotting the distribution of zeros in cells eliminating all genes without counts
#' @description This function plots the zeros distributions in cells and removes genes without counts
#' @param data.folder, a character string indicating the folder where filtered sel file will be saved. The saved file will have the prefix *filtered*.
#' @param counts.matrix, a character string indicating the the name of tab delimited file of cells un-normalized expression counts
#' @param threshold, an number from 0 to 1 indicating the fraction of max accepted zeros in each gene. 0 is set as default and it eliminates only genes which do not ave any expression in any cell.
#' @return a PDF providing zeros distributions before removal of all genes without counts. A tab delimited file with the prefix *filtered* in which the filtered data are saved.
#' @examples
#' \dontrun{
#'     #downloading fastq files
#'     system("wget http://130.192.119.59/public/singlecells_counts.txt.gz")
#'     system("gzip -d singlecells_counts.txt.gz")
#'     filterZeros(data.folder=getwd(),counts.matrix="singlecells_counts.txt", threshold=0.1)
#' }
#' @export
filterZeros <- function(data.folder=getwd(), counts.matrix, threshold=0){
  counts <- read.table(counts.matrix, sep="\t", header=T, row.names = 1)
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
      hist(counts.sum, col=rgb(1,0,0,0.5), main="",  xlab="# genes without zeros")
      hist(counts.sum0, col=rgb(0,0,1,0.5), add=T)
      box()
  dev.off()
  write.table(counts.n0, paste("filtered_",counts.matrix,sep=""), sep="\t", col.names=NA)

}

#' @title Plotting genes to UMIs relationship.
#' @description This function create a plot in which genes are plotted with respect to total counts UMI for each cell.
#' @param data.folder, a character string indicating the folder where output will be saved.
#' @param counts.matrix, a character string indicating the the name of tab delimited file  of cells un-normalized expression counts.
#' @param umiXgene, a integer defining how many UMI are required to call a gene as present. default: 3
#' @return pdf with the cells counts distributions: genes.umi.pdf
#'
#' @export
#'
#' @examples
#' \dontrun{
#'     #downloading fastq files
#'     system("wget http://130.192.119.59/public/singlecells_counts.txt.gz")
#'     system("gzip -d singlecells_counts.txt.gz")
#'     genesUmi(data.folder=getwd(), counts.matrix="singlecells_counts.txt", umiXgene=3)
#' }
genesUmi <- function(data.folder, counts.matrix, umiXgene=3){

  #running time 1
  ptm <- proc.time()
  #running time 1
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    return()
  }
  setwd(data.folder)
  tmp <- read.table(counts.matrix, sep="\t", heder=T, row.names=1, stringsAsFactors = F)
  genes <- apply(tmp,2, function(x, umiXgene){
     x[which(x <  umiXgene)] <- 0
     x[which(x >=  umiXgene)] <- 1
  })
  umi <- apply(tmp,2, sum)
  pdf("genes.umi.pdf")
     plot(log10(umi), genes, xlab="log10 UMI", ylab="# of genes")
  dev.off()
}

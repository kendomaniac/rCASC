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

  setwd(data.folder)
  tmp <- read.table(counts.matrix, sep="\t", header=T, row.names=1, stringsAsFactors = F)
  genes <- list()
  for(i in 1:dim(tmp)[2]){
    x = rep(0, dim(tmp)[1])
    x[which(tmp[,i] >=  umiXgene)] <- 1
    genes[[i]] <- x
  }
  genes <- as.data.frame(genes)
  genes.sum <-  apply(genes,2, sum)
  umi.sum <- apply(tmp,2, sum)
  pdf("genes.umi.pdf")
     plot(log10(umi.sum), genes.sum, xlab="log10 UMI", ylab="# of genes")
  dev.off()
  #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(data.folder)
  dir <- dir[grep("run.info",dir)]
  if(length(dir)>0){
    con <- file("run.info", "r")
    tmp.run <- readLines(con)
    close(con)
    tmp.run[length(tmp.run)+1] <- paste("casc checkCountDepth user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc checkCountDepth system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc checkCountDepth elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("casc user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc geneU system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc checkCountDepth elapsed run time mins ",ptm[3]/60, sep="")

    writeLines(tmp.run,"run.info")
  }

}


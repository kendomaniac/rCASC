#' @title Plotting genes to UMIs relationship.
#' @description This function create a plot in which genes are plotted with respect to total counts UMI for each cell.
#' @param file, a character string indicating the path of the file tab delimited  of cells un-normalized expression counts.
#' @param umiXgene, a integer defining how many UMI are required to call a gene as present. default: 3
#' @param sep, separator used in count file, e.g. '\\t', ','
#' @return pdf with the cells counts distributions: genes.umi.pdf
#'
#' @export
#'
#' @examples
#' \dontrun{
#'     #downloading fastq files
#'     system("wget http://130.192.119.59/public/singlecells_counts.txt.gz")
#'     system("gzip -d singlecells_counts.txt.gz")
#'     genesUmi(file=paste(getwd(),"singlecells_counts.txt",sep="/"), umiXgene=3, sep="\t")
#' }
genesUmi <- function(file, umiXgene=3, sep){

data.folder=dirname(file)
positions=length(strsplit(basename(file),"\\.")[[1]])
matrixNameC=strsplit(basename(file),"\\.")[[1]]
counts.table=paste(matrixNameC[seq(1,positions-1)],collapse="")
  matrixName=counts.table
file.type=strsplit(basename(basename(file)),"\\.")[[1]][positions]

  
  counts.matrix=matrixName

  #running time 1
  ptm <- proc.time()
  #running time 1

  setwd(data.folder)
  tmp <- read.table(paste(counts.matrix,".",file.type,sep=""), sep=sep, header=T, row.names=1, stringsAsFactors = F)
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
    tmp.run[length(tmp.run)+1] <- paste("casc geneUmi user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc geneUmi system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc geneUmi elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("casc user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc geneUmi system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc geneUmi elapsed run time mins ",ptm[3]/60, sep="")

    writeLines(tmp.run,"run.info")
  }

}


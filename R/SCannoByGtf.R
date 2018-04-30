#' @title Annotating single cell counts table using ENSEMBL gtf and refGenome CRAN package
#' @description This function executes the docker container annotate.1, where refGenome is used to annotate a single cell counts table with ensembl gene ids on first column using ENSEMBL GTF annotation
#' @param group, a character string. Two options: \code{"sudo"} or \code{"docker"}, depending to which group the user belongs
#' @param data.folder, a character string indicating where counts file is located
#' @param counts.table, a character string indicating the counts table file
#' @param gtf.name, a character string indicating the ENSEMBL gtf file
#' @param biotype, a character string the biotypes of interest
#' @param mt, a boolean to define if mitochondrial genes have to be removed, FALSE mean that mt genes are removed
#' @param ribo.proteins, a boolean to define if ribosomal proteins have to be removed, FALSE mean that ribosomal proteins (gene names starting with rpl or rps) are removed
#' @param file.type, type of file: txt tab separated columns csv comma separated columns
#' @param umiXgene,  a integer defining how many UMI are required to call a gene as present. default: 3
#' @author Raffaele Calogero

#' @return one file: annotated_counts table, where ensembl ids are linked to gene symbols

#' @import utils
#' @examples
#' \dontrun{
#'     system("wget http://130.192.119.59/public/testSCumi_mm10.csv.zip")
#'     library(casc)
#'     #filtering low quality cells
#'     lorenzFilter(group="docker",scratch.folder="/data/scratch/",
#'                  data.folder=getwd(),matrixName="filtered_testSCumi_mm10",
#'                  p_value=0.05,format="txt",separator='\t')
#'     #running annotation and removal of mit and ribo proteins genes
#'     #download mouse GTF for mm10
#'     system("wget ftp://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/Mus_musculus.GRCm38.92.gtf.gz")
#'     scannobyGtf(group="docker", data.folder=getwd(), counts.table="lorenz_filtered_testSCumi_mm10.txt",
#'                   gtf.name="Mus_musculus.GRCm38.92.gtf",
#'                   biotype="protein_coding", mt=TRUE, ribo.proteins=TRUE, file.type="txt", umiXgene=3)
#' }
#'
#' @export
scannobyGtf <- function(group=c("docker","sudo"), data.folder=getwd(), counts.table, gtf.name,
                        biotype=NULL, mt=c(TRUE, FALSE), ribo.proteins=c(TRUE, FALSE), file.type=c("txt","csv"), umiXgene=3){

  #remembering actual folder
  home <- getwd()
  #setting rsem output folder as working dir
  setwd(data.folder)

  #running time 1
  ptm <- proc.time()
  #running time 1
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    return()
  }

  if(group=="sudo"){
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",data.folder,":/data/scratch -v -d docker.io/repbioinfo/r332.2017.01 Rscript /bin/.scannoByGtf.R ", counts.table, " ", gtf.name, " ", biotype, " ", mt, " ", ribo.proteins, " ", file.type, sep="")
    resultRun <- runDocker(group="sudo",container="docker.io/repbioinfo/r332.2017.01", params=params)
  }else{
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",data.folder,":/data/scratch -v -d docker.io/repbioinfo/r332.2017.01 Rscript /bin/.scannoByGtf.R ", counts.table, " ", gtf.name, " ", biotype, " ", mt, " ", ribo.proteins, " ", file.type, sep="")
    resultRun <- runDocker(group="docker",container="docker.io/repbioinfo/r332.2017.01", params=params)
  }

  if(resultRun==0){
    cat("\nGTF based annotation is finished is finished\n")
  }
  dir <- dir(data.folder)
  files <- dir[grep(counts.table, dir)]
  files.tmp <- dir[grep(paste(file.type,'$',sep=""), dir)]
  if(length(files) == 2){
    files.annotated <- dir[grep("^annotated", dir)]
    output <- intersect(files.tmp, files.annotated)
    input <- setdiff(files.tmp, files.annotated)
    #plotting the genes vs umi all cells
    if(file.type=="txt"){
       tmp0 <- read.table(input, sep="\t", header=T, row.names=1)
       tmp <- read.table(output, sep="\t", header=T, row.names=1)
    }else{
       tmp0 <- read.table(input, sep=",", header=T, row.names=1)
       tmp <- read.table(output, sep=",", header=T, row.names=1)
    }
    genes0 <- list()
    for(i in 1:dim(tmp0)[2]){
      x = rep(0, dim(tmp0)[1])
      x[which(tmp0[,i] >=  umiXgene)] <- 1
      genes0[[i]] <- x
    }
    genes0 <- as.data.frame(genes0)
    genes.sum0 <-  apply(genes0,2, sum)
    umi.sum0 <- apply(tmp0,2, sum)

    genes <- list()
    for(i in 1:dim(tmp)[2]){
      x = rep(0, dim(tmp)[1])
      x[which(tmp[,i] >=  umiXgene)] <- 1
      genes[[i]] <- x
    }
    genes <- as.data.frame(genes)
    genes.sum <-  apply(genes,2, sum)
    umi.sum <- apply(tmp,2, sum)

    pdf("annotated_genes.pdf")
    plot(log10(umi.sum0), genes.sum0, xlab="log10 UMI", ylab="# of genes",
         xlim=c(log10(min(c(umi.sum0 + 1, umi.sum +1))), log10(max(c(umi.sum0 + 1, umi.sum + 1)))),
         ylim=c(min(c(genes.sum0, genes.sum)), max(c(genes.sum0, genes.sum))), type="n")
    points(log10(umi.sum0 + 1), genes.sum0, pch=19, cex=0.2, col="blue")
    points(log10(umi.sum + 1), genes.sum, pch=19, cex=0.2, col="red")
    legend("topleft",legend=c("All","Retained & annotated"), pch=c(15,15), col=c("blue", "red"))
    dev.off()
  }

  cat("\nannotated_genes.pdf is ready\n")

  #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(data.folder)
  dir <- dir[grep("run.info",dir)]
  if(length(dir)>0){
  con <- file("run.info", "r")
  tmp.run <- readLines(con)
  close(con)
    tmp.run[length(tmp.run)+1] <- paste("scannoByGtf user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("scannoByGtf system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("scannoByGtf elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("scannoByGtf user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("scannoByGtf system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("scannoByGtf elapsed run time mins ",ptm[3]/60, sep="")

    writeLines(tmp.run,"run.info")
  }

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", container.id, " >& ", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  system("rm -fR anno.info")
  system("rm -fR dockerID")
  system(paste("cp ",paste(path.package(package="casc"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))

  setwd(home)
}

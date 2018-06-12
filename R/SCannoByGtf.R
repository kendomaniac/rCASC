#' @title Annotating single cell counts table using ENSEMBL gtf and refGenome CRAN package
#' @description This function executes the docker container annotate.1, where refGenome is used to annotate a single cell counts table with ensembl gene ids on first column using ENSEMBL GTF annotation
#' @param group, a character string. Two options: \code{"sudo"} or \code{"docker"}, depending to which group the user belongs
#' @param file, a character string indicating the folder where input data are located and where output will be written and matrix name "/bin/users/matrix.csv"
#' @param gtf.name, a character string indicating the ENSEMBL gtf file
#' @param biotype, a character string the biotypes of interest
#' @param mt, a boolean to define if mitocondrial genes have to be removed, FALSE mean that mt genes are removed
#' @param ribo.proteins, a boolean to define if ribosomal proteins have to be removed, FALSE mean that ribosomal proteins (gene names starting with rpl or rps) are removed
#' @param umiXgene,  a integer defining how many UMI are required to call a gene as present. default: 3
#' @author Raffaele Calogero

#' @return one file: annotated_counts table, where ensembl ids are linked to gene symbols

#' @import utils
#' @examples
#' \dontrun{
#'         system("wget http://130.192.119.59/public/testSCumi_mm10.csv.zip")
#'      library(casc)
#'      system("unzip testSCumi_mm10.csv.zip")
#'      #filtering low quality cells
#'      lorenzFilter(group="docker",scratch.folder="/data/scratch/",
#'                   file=paste(getwd(),"testSCumi_mm10.csv",sep="/"),
#'                   p_value=0.05,separator=',')
#'      #running annotation and removal of mit and ribo proteins genes
#'      #download mouse GTF for mm10
#'      system("wget ftp://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/Mus_musculus.GRCm38.92.gtf.gz")
#'      system("gunzip Mus_musculus.GRCm38.92.gtf.gz")
#'      scannobyGtf(group="docker", file=paste(getwd(),"lorenz_testSCumi_mm10.csv",sep="/"),
#'                    gtf.name="Mus_musculus.GRCm38.92.gtf",
#' biotype="protein_coding", mt=TRUE, ribo.proteins=TRUE,umiXgene=3) 

#' }
#'
#' @export
scannobyGtf <- function(group=c("docker","sudo"),file, gtf.name,
                        biotype=NULL, mt=c(TRUE, FALSE), ribo.proteins=c(TRUE, FALSE), umiXgene=3){

  
  data.folder=dirname(file)
positions=length(strsplit(basename(file),"\\.")[[1]])
matrixNameC=strsplit(basename(file),"\\.")[[1]]
counts.table=paste(matrixNameC[seq(1,positions-1)],collapse="")
  matrixName=counts.table
file.type=strsplit(basename(basename(file)),"\\.")[[1]][positions]
scratch.folder=data.folder
  counts.table=paste(counts.table,".",file.type,sep="") 
 
  
  #running time 1
  ptm <- proc.time()
  #setting the data.folder as working folder
  if (!file.exists(data.folder)){
    cat(paste("\nIt seems that the ",data.folder, " folder does not exist\n"))
    return(2)
  }
  
  #storing the position of the home folder  
  home <- getwd()
  setwd(data.folder)
  #initialize status
  system("echo 0 > ExitStatusFile 2>&1")
  
  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    system("echo 10 > ExitStatusFile 2>&1") 
    setwd(home)
    return(10)
  }
  

  
  #check  if scratch folder exist
  if (!file.exists(scratch.folder)){
    cat(paste("\nIt seems that the ",scratch.folder, " folder does not exist\n"))
    system("echo 3 > ExitStatusFile 2>&1")
    setwd(data.folder)
    return(3)
  }
  tmp.folder <- gsub(":","-",gsub(" ","-",date()))
  scrat_tmp.folder=file.path(scratch.folder, tmp.folder)
  writeLines(scrat_tmp.folder,paste(data.folder,"/tempFolderID", sep=""))
  cat("\ncreating a folder in scratch folder\n")
  dir.create(file.path(scrat_tmp.folder))
  #preprocess matrix and copying files 

  #executing the docker job
 print(file.type)
 params <- paste("--cidfile ",data.folder,"/dockerID -v ",data.folder,":/data/scratch -d docker.io/repbioinfo/r332.2017.01 Rscript /bin/.scannoByGtf.R ", counts.table, " ", gtf.name, " ", biotype, " ", mt, " ", ribo.proteins, " ", file.type, sep="")
resultRun <- runDocker(group=group, params=params)
  
  #waiting for the end of the container work
  if(resultRun==0){
    cat("\nGTF based annotation is finished is finished\n")
  }
   dir <- dir(data.folder)
  files <- dir[grep(counts.table, dir)]
  if(length(files) == 2){
    files.annotated <- dir[grep("^annotated", dir)]
    output <- intersect(files, files.annotated)
    input <- setdiff(files, files.annotated)
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

    pdf(paste(matrixName,"_annotated_genes.pdf",sep=""))
    plot(log10(umi.sum0), genes.sum0, xlab="log10 UMI", ylab="# of genes",
         xlim=c(log10(min(c(umi.sum0 + 1, umi.sum +1))), log10(max(c(umi.sum0 + 1, umi.sum + 1)))),
         ylim=c(min(c(genes.sum0, genes.sum)), max(c(genes.sum0, genes.sum))), type="n")
    points(log10(umi.sum0 + 1), genes.sum0, pch=19, cex=0.2, col="blue")
    points(log10(umi.sum + 1), genes.sum, pch=19, cex=0.2, col="red")
    legend("topleft",legend=c("All","Retained & annotated"), pch=c(15,15), col=c("blue", "red"))
    
    #/////////////////////////////////////////////
    
xCoord=log10(colSums(tmp))
b=tmp
b[b<3]=0
b[b>=3]=1
yCoord=colSums(b)


xCoord2=log10(colSums(tmp0))
b2=tmp0
b2[b2<3]=0
b2[b2>=3]=1
yCoord2=colSums(b2)


 plot(yCoord2,yCoord2-yCoord,cex=0.2,pch=19,col="purple",xlab="# of genes", ylab="genesWMT&rib - genes-MT-RB")

    
    
    
    #////////////////////////////////////////////////
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
    tmp.run[length(tmp.run)+1] <- paste("user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("elapsed run time mins ",ptm[3]/60, sep="")

    writeLines(tmp.run,"run.info")
  }

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  
  

  
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
 # system(paste("cp ",paste(path.package(package="docker4seq"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}

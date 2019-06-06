#' @title A function to selectec top X on the basis of gene/transcript expression
#' @description This function select the X top genes givea a user defined threshold
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param file, a character string indicating the path of the file. IMPORTANT: full path to the file MUST be included
#' @param threshold, integer used for filtering indicate the number of top expressed genes to be selected
#' @param logged, boolean TRUE or FALSE, if FALSE gene expression data are log10 transformed before being plotted.
#' @param type, expression refers to the selection of the top expressed genes, variance to the the selectionof the top variable genes
#' @param separator, separator used in count file, e.g. '\\t', ','
#'
#' @author Raffaele Calogero, raffaele.calogero [at] unito [dot] it, UNITO
#'
#' @return a filtered tab delimited file and a histogram of the gene by gene total expression
#'
#' @examples
#'\dontrun{
#'
#'  system("wget http://130.192.119.59/public/singlecells_counts.txt.gz")
#'  system("gzip -d singlecells_counts.txt.gz")
#'  topx(group="docker", file=paste(getwd(), "singlecells_counts.txt", sep="/"),threshold=10000, logged=FALSE, type="expression", separator="\t")
#' }
#'
#' @export
topx <- function(group=c("sudo","docker"),file, threshold, separator, logged=FALSE, type=c("expression", "variance")){
  data.folder=dirname(file)
  positions=length(strsplit(basename(file),"\\.")[[1]])
  matrixNameC=strsplit(basename(file),"\\.")[[1]]
  matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
  format=strsplit(basename(basename(file)),"\\.")[[1]][positions]
  if(separator=="\t"){
    separator="tab"
  }
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
  system("echo 0 >& ExitStatusFile")
  
  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    system("echo 10 >& ExitStatusFile")
    setwd(home)
    return(10)
  }
  if(logged){logged="TRUE"}else{logged="FALSE"}
  #executing the docker job
  params <- paste("--cidfile ",data.folder,"/dockerID -v ",data.folder, ":/data -d docker.io/repbioinfo/desc.2018.01 Rscript /bin/top.R ", matrixName," ",format," ",separator, " ", logged, " ", threshold," ",type, sep="")
  resultRun <- runDocker(group=group, params=params)
  
  #waiting for the end of the container work
  if(resultRun==0){
    cat("\nData filtering is finished\n")
  }
  
  #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(data.folder)
  dir <- dir[grep("run.info",dir)]
  if(length(dir)>0){
    con <- file("run.info", "r")
    tmp.run <- readLines(con)
    close(con)
    tmp.run[length(tmp.run)+1] <- paste("topX user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("topX system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("topX elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("topX run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("topX system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("topX elapsed run time mins ",ptm[3]/60, sep="")
    
    writeLines(tmp.run,"run.info")
  }
  
  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  #removing temporary folder
  #  cat("\n\nRemoving the temporary file ....\n")
  #  system(paste("rm -fR ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}

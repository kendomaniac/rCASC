#' @title Running SCnorm  normalization
#' @description This function is a wrapper for SCnorm: robust normalization of single-cell RNA-seq data (Bacher et al. Nature Methods 2017, 14:584–586)
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param file, a character string indicating the path of the file. IMPORTANT: full path to the file MUST be included. Only tab delimited files are supported
#' @param conditions, vector of condition labels, this should correspond to the columns of the un-normalized expression matrix. If not provided data is assumed to come from same condition/batch.
#' @param outputName, specify the path and/or name of output files.
#' @param nCores, number of cores to use, default is detectCores() - 1.
#' @param filtercellNum, the number of non-zero expression estimate required to include the genes into the SCnorm fitting (default = 10). The initial grouping fits a quantile regression to each gene, making this value too low gives unstable fits.
#' @param ditherCount, FALSE of TRUE. Setting to TRUE might improve results with UMI data.
#' @param PropToUse, as default is set to 0.25, but to increase speed with large data set could be reduced, e.g. 0.1
#' @param PrintProgressPlots, producesa  plot as SCnorm determines the optimal number of groups
#' @param FilterExpression, a value indicating exclude genes having median of non-zero expression below this threshold from count-depth plots
#' @return  a tab delimited file containing the normalized data and a list of the discarded genes.
#' @examples
#' \dontrun{
#'     #downloading fastq files
#'     system("wget http://130.192.119.59/public/example_UMI.txt.zip")
#'     unzip("example_UMI.txt.zip")
#'     conditions=rep(1,12)
#'     scnorm(group="docker", file=paste(getwd(), "example_UMI.txt", sep="/"),
#'     conditions=conditions,outputName="example_UMI", nCores=8, filtercellNum=10,
#'     ditherCount=TRUE, PropToUse=0.1, PrintProgressPlots=FALSE, FilterExpression=1)

#' }
#' @export
scnorm <- function(group=c("sudo","docker"), file, conditions=NULL, outputName, nCores=8, filtercellNum = 10, ditherCount=FALSE, PropToUse=0.1, PrintProgressPlots=FALSE, FilterExpression=0){

  home <- getwd()
  data.folder=dirname(file)

  positions=length(strsplit(basename(file),"\\.")[[1]])
  matrixNameC=strsplit(basename(file),"\\.")[[1]]
  matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
  format=strsplit(basename(basename(file)),"\\.")[[1]][positions]
  if(format=="txt"){
    counts.matrix <- paste(matrixName, format, sep=".")
  }else{
    cat("\nOnly tab delimited files with extention txt are supported\n")
    return("Not a tab delimited file")
  }

  #running time 1
  ptm <- proc.time()
  #running time 1
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    return()
  }

  setwd(data.folder)
  if(is.null(conditions)){
    cat("\nERROR: Conditions are missing\n")
  }else{
    conditions <- paste(conditions, collapse = "_")
  }

  #checking eligibility to scnorm
  tmp <- read.table(paste(data.folder, counts.matrix, sep="/"), sep="\t", header=T, row.names=1)
  tmp.sum <- apply(tmp,2,sum)
  if(min(tmp.sum) < 10000){

    cat("\n There are ",length(which(tmp.sum) < 10000), " samples out of ", length(tmp.sum), " with less than 10000 counts\n")
    cat("Remove them before applying scnorm\n")
    return("\nERROR: some of the cells have less than 10000 counts\n")

  }

  params <- paste("--cidfile ",data.folder,"/dockerID -v ", data.folder,":/data -d docker.io/repbioinfo/scnorm.2018.01 Rscript /bin/scnorm.R ",counts.matrix," ",conditions," ",outputName," ",nCores," ",filtercellNum, " ",ditherCount," ",PropToUse," ", PrintProgressPlots," ", FilterExpression, sep="")
  resultRun <- runDocker(group=group, params=params)

  if(resultRun==0){
    cat("\nSCnorm is finished\n")
  }

  #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(data.folder)
  dir <- dir[grep("run.info",dir)]
  if(length(dir)>0){
    con <- file("run.info", "r")
    tmp.run <- readLines(con)
    close(con)
    tmp.run[length(tmp.run)+1] <- paste("casc SCnorm user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc SCnorm system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc SCnorm elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("casc user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc SCnorm system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc SCnorm elapsed run time mins ",ptm[3]/60, sep="")

    writeLines(tmp.run,"run.info")
  }

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", substr(container.id,1,12), " >& ", "SCnorm_", substr(container.id,1,12),".log", sep=""))

  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system("rm -fR anno.info")
  system("rm -fR dockerID")
  system(paste("cp ",paste(path.package(package="casc"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))

  system(paste("docker rm ", container.id, sep=""))
  setwd(home)

}

#' @title Running SCnorm  checkCountDepth test.
#' @description This function executes a check on the data count-depth relationship used by SCnorm.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs.
#' @param file, a character string indicating the path of the file. IMPORTANT: full path to the file MUST be included. Only tab delimited files are supported
#' @param conditions, vector of condition labels, this should correspond to the columns of the un-normalized expression matrix. If not provided data is assumed to come from same condition/batch.
#' @param FilterCellProportion, a value indicating the proportion of non-zero expression estimates required to include the genes into the evaluation. Default is .10, and will not go below a proportion which uses less than 10 total cells/samples
#' @param FilterExpression, a value indicating exclude genes having median of non-zero expression below this threshold from count-depth plots
#' @param ditherCounts, Setting to TRUE might improve results with UMI data, default is FALSE
#' @param outputName, specify the path and/or name of output files.
#' @param nCores, number of cores to use, default is detectCores() - 1.
#' @return pdf with the cells counts distributions
#' @examples
#' \dontrun{
#'     #UMI 3' end analysis
#'     system("wget http://130.192.119.59/public/example_UMI.txt.zip")
#'     unzip("example_UMI.txt.zip")
#'     conditions=rep(1,12)
#'     checkCountDepth(group="docker", file=paste(getwd(), "example_UMI.txt", sep="/"),
#'     conditions=conditions, FilterCellProportion=0.1, FilterExpression=0,
#'     ditherCounts=TRUE, outputName="example_UMI", nCores=8)
#'
#' }
#' @export
checkCountDepth <- function(group=c("sudo","docker"), file, conditions=NULL, FilterCellProportion=0.1, FilterExpression=0, ditherCounts=FALSE, outputName, nCores=8){

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

  home <- getwd()
  setwd(data.folder)
  if(is.null(conditions)){
    cat("\nERROR: Conditions are missing\n")
  }else{
    conditions <- paste(conditions, collapse = "_")
  }

  params <- paste("--cidfile ",data.folder,"/dockerID -v ", data.folder,":/data -d docker.io/repbioinfo/scnorm.2018.01 Rscript /bin/checkCountDepth.R ",counts.matrix," ",conditions," ", FilterCellProportion, " ", FilterExpression," ", ditherCounts, " ", outputName," ",nCores, sep="")
  resultRun <- runDocker(group=group, params=params)

#  if(group=="sudo"){
#    params <- paste("--cidfile ",data.folder,"/dockerID -v ", data.folder,":/data -d docker.io/repbioinfo/scnorm.2018.01 Rscript /bin/checkCountDepth.R ",counts.matrix," ",conditions," ", FilterCellProportion, " ", FilterExpression," ", ditherCounts, " ", outputName," ",nCores, sep="")
#    resultRun <- runDocker(group="sudo",container="docker.io/repbioinfo/scnorm.2018.01", params=params)
#  }else{
#    params <- paste("--cidfile ",data.folder,"/dockerID -v ", data.folder,":/data -d docker.io/repbioinfo/scnorm.2018.01 Rscript /bin/checkCountDepth.R ",counts.matrix," ",conditions," ",FilterCellProportion, " ", FilterExpression," ", ditherCounts, " ", outputName," ",nCores, sep="")
#    resultRun <- runDocker(group="docker",container="docker.io/repbioinfo/scnorm.2018.01", params=params)
#  }

  if(resultRun==0){
    cat("\ncheckCountDepth test is finished\n")
  }

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
    tmp.run[length(tmp.run)+1] <- paste("casc checkCountDepth system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc checkCountDepth elapsed run time mins ",ptm[3]/60, sep="")

    writeLines(tmp.run,"run.info")
  }

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  #command below makes an error!
  system(paste("docker logs ", substr(container.id,1,12), " >& ", "checkCountDepth_", substr(container.id,1,12),".log", sep=""))
  # system(paste("docker rm ", container.id, sep=""))

  #removing temporary folder
  cat("\n\nRemoving the checkCountDepth temporary file ....\n")
  system(paste("rm  -f ",data.folder,"/dockerID", sep=""))
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)

}

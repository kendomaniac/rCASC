#' @title Running SCnorm  normalization
#' @description This function is a wrapper for SCnorm: robust normalization of single-cell RNA-seq data (Bacher et al. Nature Methods 2017, 14:584–586)
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param file, a character string indicating the path of the file. IMPORTANT: full path to the file MUST be included. Only tab delimited files are supported
#' @param outputName, specify the path and/or name of output files.
#' @param normMethod, a string identifying the normalization method: CLR_FN, DESEQ_FN, FQ_FN, SCRAN_FN, SUM_FN, TMM_FN, UQ_FN, more info on the vignette.
#' @return a tab delimited file containing the normalized data.
#' @examples
#' \dontrun{
#'     #downloading fastq files
#'     system("wget http://130.192.119.59/public/example_UMI.txt.zip")
#'     unzip("example_UMI.txt.zip")
#'     umiNorm(group="docker", file=paste(getwd(), "example_UMI.txt", sep="/"),
#'     outputName="example_UMI", normMethod="TMM_FN")
#' }
#' @export
umiNorm <- function(group=c("sudo","docker"), file, outputName, normMethod=c("CLR_FN", "DESEQ_FN", "FQ_FN", "SCRAN_FN", "SUM_FN", "TMM_FN", "UQ_FN")){

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

  params <- paste("--cidfile ",data.folder,"/dockerID -v ", data.folder,":/data -d docker.io/repbioinfo/scone.2018.01 Rscript /bin/scone.R ",counts.matrix," ", paste(outputName, format, sep=".")," ", normMethod, sep="")
  resultRun <- runDocker(group=group, params=params)

#  if(group=="sudo"){
#    params <- paste("--cidfile ",data.folder,"/dockerID -v ", data.folder,":/data -d docker.io/repbioinfo/scone.2018.01 Rscript /bin/scone.R ",counts.matrix," ", outputName," ", normMethod, sep="")
#    resultRun <- runDocker(group="sudo",container="docker.io/repbioinfo/scone.2018.01", params=params)
#  }else{
#    params <- paste("--cidfile ",data.folder,"/dockerID -v ", data.folder,":/data -d docker.io/repbioinfo/scone.2018.01 Rscript /bin/scone.R ",counts.matrix," ", outputName," ",normMethod, sep="")
#    resultRun <- runDocker(group="docker",container="docker.io/repbioinfo/scone.2018.01", params=params)
#  }

  if(resultRun==0){
    cat("\nScone is finished\n")
  }

  #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(data.folder)
  dir <- dir[grep("run.info",dir)]
  if(length(dir)>0){
    con <- file("run.info", "r")
    tmp.run <- readLines(con)
    close(con)
    tmp.run[length(tmp.run)+1] <- paste("casc scone user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc scone system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc scone elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("casc user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc scone system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc scone elapsed run time mins ",ptm[3]/60, sep="")

    writeLines(tmp.run,"run.info")
  }

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", substr(container.id,1,12), " >& ", "scone_", substr(container.id,1,12),".log", sep=""))

  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system("rm -fR anno.info")
  system("rm -fR dockerID")
  system(paste("cp ",paste(path.package(package="casc"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))

  system(paste("docker rm ", container.id, sep=""))
  setwd(home)
}

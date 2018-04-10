#' @title Running SCnorm  normalization
#' @description This function is a wrapper for SCnorm: robust normalization of single-cell RNA-seq data (Bacher et al. Nature Methods 2017, 14:584–586)
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param data.folder, a character string indicating the folder where comma separated file of cells log10 counts is saved
#' @param counts.matrix, a character string indicating the the name of tab delimited file  file of cells un-normalized expression counts
#' @param outputName, specify the path and/or name of output files.
#' @param normMethod, a string identifying the normalization method: CLR_FN, DESEQ_FN, FQ_FN, SCRAN_FN, SUM_FN, TMM_FN, UQ_FN
#' @return a tab delimited file containing the normalized data.
#' @examples
#' \dontrun{
#'     #downloading fastq files
#'     system("wget http://130.192.119.59/public/example_UMI.txt.zip")
#'     unzip("example_UMI.txt.zip")
#'     conditions=rep(1,12)
#'     umiNorm(group="docker", data.folder=getwd(), counts.matrix="example_UMI.txt",
#'     outputName="example_UMI", normMethod="TMM_FN")
#' }
#' @export
umiNorm <- function(group=c("sudo","docker"), data.folder=getwd(), counts.matrix, outputName, normMethod=c("CLR_FN", "DESEQ_FN", "FQ_FN", "SCRAN_FN", "SUM_FN", "TMM_FN", "UQ_FN")){

  #running time 1
  ptm <- proc.time()
  #running time 1
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    return()
  }
  setwd(data.folder)

  if(group=="sudo"){
    params <- paste("--cidfile ",data.folder,"/dockerID -v ", data.folder,":/data -d docker.io/repbioinfo/scone.2018.01 Rscript /bin/scnorm.R ",counts.matrix," ", outputName," ", normMethod, sep="")
    resultRun <- runDocker(group="sudo",container="docker.io/repbioinfo/scone.2018.01", params=params)
  }else{
    params <- paste("--cidfile ",data.folder,"/dockerID -v ", data.folder,":/data -d docker.io/repbioinfo/scone.2018.01 Rscript /bin/scnorm.R ",counts.matrix," ", outputName," ",normMethod, sep="")
    resultRun <- runDocker(group="docker",container="docker.io/repbioinfo/scone.2018.01", params=params)
  }

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

}

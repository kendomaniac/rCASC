#' @title Imputing dropouts
#' @description This function executes a wrapper for scImpute: Accurate And Robust Inputation For Single Cell RNA-Seq Data (Li WV Nature Communications, vol. 9, Article number: 997, 2018)
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param data.folder, a character string indicating the folder where tab separated file of cells counts is located
#' @param counts.matrix, a character string indicating the  tab separated file of cells  counts
#' @param drop.thre, A number between 0 and 1, specifying the threshold to determine dropout values
#' @param cores, a integer specifying the number of cores used for parallel computation.
#' @param refining, boolean FALSE, TRUE, if TRUE it execute again scImpute but setting a different drop.thre. Rerunning to change drop.thre hase to be done in the same folder where the first run was done.
#' @return A matrix file with imputed data.
#' @examples
#' \dontrun{
#'     #downloading fastq files
#'     system("wget http://130.192.119.59/public/singlecells_counts.txt.gz")
#'     system("gzip -d singlecells_counts.txt.gz")
#'     cascImpute(group="docker", data.folder=getwd(),
#'          counts.matrix="singlecells_counts.txt",
#'          drop.thre=0.5, cores=8, refining=FALSE)
#' }
#'
#' @export
impute <- function(group=c("sudo","docker"), data.folder, counts.matrix, drop.thre, cores, refining=FALSE){

  if(group=="sudo"){
    params <- paste("--cidfile ",data.folder,"/dockerID -v ", data.folder, ":/data -d docker.io/repbioinfo/r340.2017.01 Rscript /bin/scimpute.R ",counts.matrix," ",drop.thre," ",cores," ", refining, sep="")
    runDocker(group="sudo",container="docker.io/repbioinfo/r340.2017.01", params=params)
  }else{
    params <- paste("--cidfile ",data.folder,"/dockerID -v ", data.folder, ":/data -d docker.io/repbioinfo/r340.2017.01 Rscript /bin/scimpute.R ",counts.matrix," ",drop.thre," ",cores," ", refining, sep="")
    runDocker(group="docker",container="docker.io/repbioinfo/r340.2017.01", params=params)
  }

  if(resultRun=="false"){
    cat("\nscImpute is finished\n")
  }

  #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(data.folder)
  dir <- dir[grep("run.info",dir)]
  if(length(dir)>0){
    con <- file("run.info", "r")
    tmp.run <- readLines(con)
    close(con)
    tmp.run[length(tmp.run)+1] <- paste("casc impute user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc impute system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc impute elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("casc impute user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc impute system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("casc impute elapsed run time mins ",ptm[3]/60, sep="")

    writeLines(tmp.run,"run.info")
  }

  #converting table
  tmp <- read.table("scimpute_count.txt", sep = " ", stringsAsFactors = FALSE ,header = T, check.names = FALSE, row.names = 1)
  write.table(tmp, "scimpute_count_tab.txt", sep="\t", col.names=NA, quote = FALSE)

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", container.id, " >& ", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system("rm -fR anno.info")
  system("rm -fR dockerID")
  system(paste("cp ",paste(path.package(package="docker4seq"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))


}

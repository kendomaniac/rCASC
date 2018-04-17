#' @title Annotating single cell counts table using ENSEMBL gtf and refGenome CRAN package
#' @description This function executes the docker container annotate.1, where refGenome is used to annotate a single cell counts table with ensembl gene ids on first column using ENSEMBL GTF annotation
#' @param group, a character string. Two options: \code{"sudo"} or \code{"docker"}, depending to which group the user belongs
#' @param data.folder, a character string indicating where counts file is located
#' @param counts.table, a character string indicating the counts table file
#' @param gtf.name, a character string indicating the ENSEMBL gtf file
#' @param biotype, a character string the biotypes of interest
#' @param mt, a boolean to define if mitocondrial genes have to be removed
#' @param file.type, type of file: txt tab separated columns csv comma separated columns
#' @author Raffaele Calogero

#' @return one file: annotated_counts table, where ensembl ids are linked to gene symbols

#' @import utils
#' @examples
#' \dontrun{
#'     #downloading fastq files
#'     system("wget http://130.192.119.59/public/genes.results.gz")
#'     system("gzip -d genes.results.gz")
#'     #running rsemannoByGtf
#'     SCannoByGtf(group="docker", data.folder=getwd(), counts.table="GSM2833284_Naive_WT_Rep1.csv",
#'                   gtf.name="Mus_musculus.GRCm38.92.gtf",
#'                   biotype, mt=TRUE, file.type="csv")
#' }
#'
#' @export
SCannoByGtf <- function(group=c("docker","sudo"), data.folder=getwd(), counts.table, gtf.name,
                        biotype, mt=c(TRUE, FALSE), file.type=c("txt","csv")){

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
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",data.folder,":/data/scratch -v -d docker.io/repbioinfo/r332.2017.01 Rscript /bin/.scannoByGtf.R ", counts.table, " ", gtf.name, " ", biotype, " ", mt, " ", file.type, sep="")
    resultRun <- runDocker(group="sudo",container="docker.io/repbioinfo/r332.2017.01", params=params)
  }else{
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",data.folder,":/data/scratch -v -d docker.io/repbioinfo/r332.2017.01 Rscript /bin/.scannoByGtf.R ", counts.table, " ", gtf.name, " ", biotype, " ", mt, " ", file.type, sep="")
    resultRun <- runDocker(group="docker",container="docker.io/repbioinfo/r332.2017.01", params=params)
  }

  if(resultRun==0){
    cat("\nGTF based annotation is finished is finished\n")
  }

  #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(rsem.folder)
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
  system(paste("cp ",paste(path.package(package="docker4seq"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))

  setwd(home)
}

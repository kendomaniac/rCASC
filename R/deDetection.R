#' @title A function allowing the identification of differentially expressed genes.
#' @description This function executes in a docker edgeR for the idnetification of differentially expressed genes in single-cells RNAseq
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param data.folder, a character string indicating the folder where input data are located and where output will be written
#' @param counts.table, a character string indicating the counts table file. IMPORTANT in the header of the file the covariate group MUST be associated to the column name using underscore, e.g. cell1_cov1
#' @param file.type, type of file: txt tab separated columns csv comma separated columns
#' @author Raffaele Calogero, raffaele.calogero [at] unito [dot] it, University of Torino, Italy
#'
#' @examples
#' \dontrun{
#'     #running deDetection
#'     system("wget http://130.192.119.59/public/buettner_counts_noSymb.txt.zip")
#'     unzip("buettner_counts_noSymb.txt.zip")
#'     lorenzFilter(group="docker", scratch.folder="/data/scratch/",
#'                 data.folder=getwd(), matrixName="buettner_counts_noSymb",
#'                 p_value=0.05, format="txt", separator='\t')
#'
#'     system("wget ftp://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/Mus_musculus.GRCm38.92.gtf.gz")
#'     system("gzip -d Mus_musculus.GRCm38.92.gtf.gz")
#'     scannobyGtf(group="docker", data.folder=getwd(),
#'                  counts.table="lorenz_buettner_counts_noSymb.txt",
#'                  gtf.name="Mus_musculus.GRCm38.92.gtf",
#'                  biotype="protein_coding", mt=FALSE, ribo.proteins=FALSE,
#'                  file.type="txt", umiXgene=3)
#'
#'     deDetection(group="docker", data.folder=getwd(),
#'                counts.table="annotated_lorenz_buettner_counts_noSymb.txt",
#'                file.type="txt")
#' }
#'
#' @export
deDetection <- function(group=c("sudo","docker"), data.folder, counts.table, file.type=c("txt","csv")){




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

  #executing the docker job
  params <- paste("--cidfile ",data.folder,"/dockerID -v ", data.folder, ":/data -d docker.io/repbioinfo/desc.2018.01 Rscript /bin/desc.R ", counts.table, " ", file.type, sep="")
  resultRun <- runDocker(group=group, params=params)

  #waiting for the end of the container work
  if(resultRun==0){
    cat("\nDifferential expression analysis is finished\n")
  }
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
    tmp.run[1] <- paste("DE detection run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("DE detection system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("DE detection elapsed run time mins ",ptm[3]/60, sep="")

    writeLines(tmp.run,"run.info")
  }

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system("rm -fR dockerID")
  system(paste("cp ",paste(path.package(package="casc"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}

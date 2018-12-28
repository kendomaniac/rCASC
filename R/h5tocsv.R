#' @title h5 to csv
#' @description This function takes h5 file from cellranger output and convert it in csv table.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param file,  path of the h5 file. Full path MUST be included.
#' 
#' @author Raffaele Calogero, raffaele [dot] calogero [at] unito [dot] com, University of Torino
#'
#' @return a csv file results_cellranger
#'
#' @examples
#' \dontrun{
#' home <- getwd()
#' library(rCASC)
#' #download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106268 the GSE106268_RAW.tar file
#' system("tar xvf GSE106268_RAW.tar")
#' 
#' h5tocvs(group="docker", file=paste(getwd(),"GSM2833284_Naive_WT_Rep1.h5",sep="/"))
#' }
#'
#'
#' @export
 
h5tocvs <- function(group=c("sudo","docker"),  file){

  dockerImage="docker.io/repbioinfo/cellranger.2018.03"
  
  data.folder=dirname(file)
  positions=length(strsplit(basename(file),"\\.")[[1]])
  matrixNameC=strsplit(basename(file),"\\.")[[1]]
  counts.table=paste(matrixNameC[seq(1,positions-1)],collapse="")
  matrixName=paste(counts.table, ".csv", sep="")

  home <- getwd()
  #running time 1
  ptm <- proc.time()
  
  #setting the data.folder as working folder
  if (!file.exists(data.folder)){
    cat(paste("\nIt seems that the ",data.folder, " folder does not exist\n"))
    return(2)
  }
  setwd(data.folder)
  
  #initialize status
  system("echo 0 > ExitStatusFile 2>&1")
  
  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    system("echo 0 > ExitStatusFile 2>&1")
    setwd(home)
    return(10)
  }
  

  params <- paste("--cidfile ", data.folder,"/dockerID -v ", data.folder, ":/data -d ", dockerImage, " /bin/cellranger mat2csv /data/ ",matrixName, sep="")
  
  #Run docker
  resultRun <- runDocker(group=group, params=params)
  #waiting for the end of the container work
  if(resultRun==0){
    system(paste("sed \'s|,|\t|g\' ", data.folder,"/",counts.table,".csv > ", data.folder,"/",counts.table,".txt", sep=""))
    cat("\nCellranger h5 matrix conversion is finished\n")
  }
  
  #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(data.folder)
  dir <- dir[grep("run.info",dir)]
  if(length(dir)>0){
    con <- file("run.info", "r")
    tmp.run <- readLines(con)
    close(con)
    tmp.run[length(tmp.run)+1] <- paste("cellranger mat2csv user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("cellranger mat2csv system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("cellranger mat2csv elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("cellranger mat2csv run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("cellranger mat2csv system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("cellranger mat2csv elapsed run time mins ",ptm[3]/60, sep="")
    
    writeLines(tmp.run,"run.info")
  }
  
  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/mat2csv_", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  system("rm -fR dockerID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}

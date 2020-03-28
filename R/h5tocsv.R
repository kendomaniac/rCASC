#' @title sparse to dense
#' @description This function takes h5 or mtx file from cellranger output and convert it in a csv table.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param file,  path of the sparse matrix file. For h5 file the full path MUST be included. For mtx matrix the folder MUST contain tsv and mtx files and the FULL path to mtx matrix MUST be provided
#' @param type, h5 refers to h5 files and 10xgenomics to the folder containing barcodes.tsv, genes.tsv and matrix.mtx
#' 
#' @author Raffaele Calogero, raffaele [dot] calogero [at] unito [dot] com, University of Torino
#'
#' @return a dense matrix, with 0s, in csv and txt format
#'
#' @examples
#' \dontrun{
#' home <- getwd()
#' library(rCASC)
#' #download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106268 the GSE106268_RAW.tar file
#' system("tar xvf GSE106268_RAW.tar")
#' h5tocsv(group="docker", file=paste(getwd(),"GSM2833284_Naive_WT_Rep1.h5",sep="/"), type="h5"))
#' system("wget http://130.192.119.59/public/annotated_setPace_10000_noC5.txt.zip")
#' unzip("annotated_setPace_10000_noC5.txt.zip")
#' csvToSparse(group="docker", scratch="/data/scratch", file=paste(getwd(), 
#'             "annotated_setPace_10000_noC5.txt", sep="/"), separator="\t")
#'             
#' h5tocsv(group="docker", file=paste(getwd(),"matrix.mtx",sep="/"), type="10xgenomics")
#' }
#'
#'
#' @export
 
h5tocsv <- function(group=c("sudo","docker"),  file, type=c("h5","10xgenomics")){

  dockerImage="docker.io/repbioinfo/cellranger.2018.03"
  
  if(type=="h5"){
    data.folder=dirname(file)
    positions=length(strsplit(basename(file),"\\.")[[1]])
    matrixNameC=strsplit(basename(file),"\\.")[[1]]
    counts.table=paste(matrixNameC[seq(1,positions-1)],collapse="")
    matrixName=paste(counts.table, ".h5", sep="")
    matrixNameOut=paste(counts.table, ".csv", sep="")
  }else{
    data.folder=dirname(file)
    matrixName=basename(file)
	matrixName = sapply(strsplit(matrixName, '\\.'), function(x)x[1])
    matrixNameOut=paste(matrixName, ".csv", sep="")
  }

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
  
  if(type =="h5"){
    params <- paste("--cidfile ", data.folder,"/dockerID -v ", data.folder, ":/data -d ", dockerImage, " /bin/cellranger mat2csv /data/",matrixName, " /data/",matrixNameOut, sep="")
  }else{
    params <- paste("--cidfile ", data.folder,"/dockerID -v ", data.folder, ":/data -d ", dockerImage, " /bin/cellranger mat2csv /data/", " /data/",matrixNameOut, sep="")
  }
  
  #Run docker
  resultRun <- runDocker(group=group, params=params)
  #waiting for the end of the container work
  if(resultRun==0 && type=="h5"){
  #  system(paste("sed \'s|,|\t|g\' ", data.folder,"/",counts.table,".csv > ", data.folder,"/",counts.table,".txt", sep=""))
    cat("\nCellranger h5 matrix conversion is finished\n")
  }else{
  #  system(paste("sed \'s|,|\t|g\' ", data.folder,"/",matrixName,".csv > ", data.folder,"/",matrixName,".txt", sep=""))
    cat("\nCellranger tsv, mtx matrices conversion is finished\n")
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

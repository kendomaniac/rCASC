 
#' @title bulkClusters
#' @description The present function create pseudo-bulk matrix from clustering.output file. The output are three files: _bulklog2, which is not normalized, _bulkColumn, which is z-scoere calculated over each column, _bulkRow, which is z-score calculated over each row
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, a character string indicating the path of the file, with file name and extension included
#' @param cl, name and path of the file clustering.output previously generated from clustering algorithm from rcasc
#' @param separator, matrix separator, ',', '\\t'
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @examples
#' \dontrun{
#'  bulkClusters(group="docker", scratch.folder="/home/user/scratch", 
#'  file="/home/user/temp/setA.csv",separator=",",
#'  cl="/home/user/temp/Results/setA/3/setA_clustering.output.csv")
#'  
#'}
#' @export
bulkClusters <- function(group=c("sudo","docker"), scratch.folder, file,separator,cl){

data.folder=dirname(file)
positions=length(strsplit(basename(file),"\\.")[[1]])
matrixNameC=strsplit(basename(file),"\\.")[[1]]
matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
format=strsplit(basename(basename(file)),"\\.")[[1]][positions]




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
  system("echo 0 > ExitStatusFile 2>&1")

  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    system("echo 10 > ExitStatusFile 2>&1")
    setwd(home)
    return(10)
  }



  #check  if scratch folder exist
  if (!file.exists(scratch.folder)){
    cat(paste("\nIt seems that the ",scratch.folder, " folder does not exist\n"))
    system("echo 3 > ExitStatusFile 2>&1")
    setwd(data.folder)
    return(3)
  }
  tmp.folder <- gsub(":","-",gsub(" ","-",date()))
  scrat_tmp.folder=file.path(scratch.folder, tmp.folder)
  writeLines(scrat_tmp.folder,paste(data.folder,"/tempFolderID", sep=""))
  cat("\ncreating a folder in scratch folder\n")
  dir.create(file.path(scrat_tmp.folder))
  #preprocess matrix and copying files

if(separator=="\t"){
separator="tab"
}

system(paste("cp ",data.folder,"/",matrixName,".",format," ",scrat_tmp.folder,"/",sep=""))
system(paste("cp ",cl," ",scrat_tmp.folder,"/",sep=""))




  #executing the docker job
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ", data.folder, ":/data -d repbioinfo/mergematrix Rscript /home/main4.R ",matrixName," ",format," ",separator," ",basename(cl),sep="")
    print(params)
resultRun <- runDocker(group=group, params=params)

  #waiting for the end of the container work
  if(resultRun==0){

    #system(paste("cp ", scrat_tmp.folder, "/* ", data.folder, sep=""))
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
    tmp.run[1] <- paste("run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("elapsed run time mins ",ptm[3]/60, sep="")

    writeLines(tmp.run,"run.info")
  }

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker rm ", container.id, sep=""))


  #Copy result folder
  cat("Copying Result Folder")
  system(paste("cp -r ",scrat_tmp.folder,"/* ",data.folder,sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}

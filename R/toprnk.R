#' @title toprnk
#' @description This function execute toprnk analysis which search for corrispondence betwee clusters of two different experiments using pseudobulk zscored on rows and the cluster specific genes from comet analysis. Thsu the function clustersBulk and cometsc have to be run in the two dataset to be integrated.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param fileX, a character string indicating the path of the pseudobulkRow file, with file name and extension included. 
#' @param fileY, a character string indicating the path of the pseudobulkRow file, with file name and extension included. 
#' @param separatorX, separator used in count file, e.g. '\\t', ','
#' @param separatorY, separator used in count file, e.g. '\\t', ','
#' @param xCometFolder, path of Comet results from X experiment
#' @param yCometFolder, path of Comet results from Y experiment
#' @param threshold, Pearson threshold
#' @param top.ranked, number of top comet genes to be used
#' @param outputFolder, where results are placed
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return 
#' @examples
#' \dontrun{
#'  

#'}
#' @export
toprnk <- function(group=c("sudo","docker"), scratch.folder, fileX,fileY, separatorX,separatorY,xCometFolder,yCometFolder,threshold=0.8,top.ranked=80,outputFolder){


data.folder=outputFolder
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

if(separatorX=="\t"){
separatorX="tab"
}

if(separatorY=="\t"){
separatorY="tab"
}

system(paste("cp ",fileX," ",scrat_tmp.folder,"/",sep=""))
fileX.file1 <- sapply(strsplit(basename(fileX), "\\."), function(x)x[1])
fileX.file2 <- sapply(strsplit(basename(fileX), "\\."), function(x)x[2])
fileX.file <- paste("X", fileX.file2, sep=".")
system(paste("mv ",scrat_tmp.folder,"/", basename(fileX)," ",scrat_tmp.folder,"/",fileX.file, sep=""))

system(paste("cp ",fileY," ",scrat_tmp.folder,"/",sep=""))
fileY.file1 <- sapply(strsplit(basename(fileY), "\\."), function(x)x[1])
fileY.file2 <- sapply(strsplit(basename(fileY), "\\."), function(x)x[2])
fileY.file <- paste("Y", fileY.file2, sep=".")
system(paste("mv ",scrat_tmp.folder,"/", basename(fileY)," ",scrat_tmp.folder,"/",fileY.file, sep=""))

system(paste("cp -r ",xCometFolder," ",scrat_tmp.folder,"/",sep=""))
system(paste("mv ", scrat_tmp.folder,"/outputdata ", scrat_tmp.folder,"/Xoutputdata", sep=""))
system(paste("cp -r ",yCometFolder," ",scrat_tmp.folder,"/",sep=""))
system(paste("mv ", scrat_tmp.folder,"/outputdata ", scrat_tmp.folder,"/Youtputdata", sep=""))

fileX=paste("/scratch/",fileX.file,sep="")
fileY=paste("/scratch/",fileY.file,sep="")
xCometFolder="/scratch/Xoutputdata"
yCometFolder="/scratch/Youtputdata"

  #executing the docker job


    params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch:Z -v ", data.folder, ":/data -d docker.io/repbioinfo/combinetoprnk Rscript /home/combineNT1NT2_toprnk.R ",top.ranked," ",fileX," ",xCometFolder," ",fileY," ",yCometFolder," ",threshold," ",separatorX," ",separatorY,sep="")
resultRun <- runDocker(group=group, params=params)

  #waiting for the end of the container work
  if(resultRun==0){
    #system(paste("cp ", scrat_tmp.folder, "/*top* ", data.folder, sep=""))
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
  system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))


  #Copy result folder
  cat("Copying Result Folder")
  system(paste("cp -r ",scrat_tmp.folder,"/*top* ",data.folder,"/",sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}
 

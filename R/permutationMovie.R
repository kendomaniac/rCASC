#' @title Permutation Movie
#' @description This function executes a ubuntu docker that create a video showing the cluster changing of all the cells
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param data.folder, a character string indicating the folder where input data are located and where output will be written
#' @param matrixName, counts table name. Matrix data file must be in data.folder. The file MUST contain RAW counts, without any modification, such as log transformation, normalizatio etc.
#' @param nCluster, number of Cluster used in Kmeans to generate the clusters that you want to merge
#' @param format, count matrix format "csv", "txt"..
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param framePP, Number of frame for each permutation
#' @param permutationNumber, Number of random permutation,have to be less or the same then the total permutation
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return Csv file with correct cluster name
#' @examples
#' \dontrun{
#'  system("wget http://130.192.119.59/public/permutationmovie.zip")
#'  unzip("permutationmovie.zip")
#'  setwd("./permutationmovie")
#'  library("CASC")
#'  permutationMovie("docker",scratch.folder="/data/scratch", data.folder=getwd(),
#'         matrixName="lorenz_Buettner", nCluster=4, format="csv", separator=",",
#'         framePP=200, permutationNumber=3)
#'}
#' @export
permutationMovie <- function(group=c("sudo","docker"), scratch.folder, data.folder,matrixName,nCluster,format,separator,framePP,permutationNumber){



  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    return()
  }
  #storing the position of the home folder
  home <- getwd()

  #running time 1
  ptm <- proc.time()
  #setting the data.folder as working folder
  if (!file.exists(data.folder)){
    cat(paste("\nIt seems that the ",data.folder, " folder does not exist\n"))
    return(2)
  }
  setwd(data.folder)
  #check  if scratch folder exist
  if (!file.exists(scratch.folder)){
    cat(paste("\nIt seems that the ",scratch.folder, " folder does not exist\n"))
    return(3)
  }
  tmp.folder <- gsub(":","-",gsub(" ","-",date()))
  scrat_tmp.folder=file.path(scratch.folder, tmp.folder)
  writeLines(scrat_tmp.folder,paste(data.folder,"/tempFolderID", sep=""))
  cat("\ncreating a folder in scratch folder\n")
  dir.create(file.path(scrat_tmp.folder))

if(separator=="\t"){
separator="tab"
}

system(paste("cp -r ",data.folder,"/Results/* ",scrat_tmp.folder,sep=""))


  #executing the docker job
  if(group=="sudo"){
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ", data.folder, ":/data -d docker.io/rcaloger/permutationmovie Rscript /home/main.R ",matrixName," ",nCluster," ",format," ",separator," ",framePP," ",permutationNumber,sep="")
    resultRun <- runDocker(group="sudo",container="docker.io/rcaloger/permutationmovie", params=params)
  }else{
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ", data.folder, ":/data -d docker.io/rcaloger/permutationmovie Rscript /home/main.R ",matrixName," ",nCluster," ",format," ",separator," ",framePP," ",permutationNumber,sep="")
    resultRun <- runDocker(group="docker",container="docker.io/rcaloger/permutationmovie", params=params)
  }
  #waiting for the end of the container work
  if(resultRun=="false"){
    #system(paste("cp -r ", scrat_tmp.folder, "/* ", data.folder,"Results", sep=""))
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
  #removing temporary folder
  cat("Copying Result Folder")
  system(paste("rm -r ",scrat_tmp.folder,"/graph ",sep=""))
  system(paste("cp -r ",scrat_tmp.folder,"/* ",data.folder,"/Results",sep=""))
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="docker4seq"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}

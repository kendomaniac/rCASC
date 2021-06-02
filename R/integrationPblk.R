#' @title integrationPsblk
#' @description This function execute integrationPsblk analysis which search for correspondence between clusters of two different experiments using clusters-pseudobulks, zscored on rows, and a subset of randomly selexte genes. Thus, the function clustersBulk has to be run in the two datasets before their comparison.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param fileX, a character string indicating the path of the pseudobulkRow file, with file name and extension included. 
#' @param fileY, a character string indicating the path of the pseudobulkRow file, with file name and extension included. 
#' @param separatorX, separator used in count file, e.g. '\\t', ','
#' @param separatorY, separator used in count file, e.g. '\\t', ','
#' @param max.genes, MAX number of random genes to be used for each cluster, default 500
#' @param outputFolder, where results are placed
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return A folder called XYpb with all the results generated. the final frequence table is saved in final_score.csv. The function produces an integrated output combining 7 thresholds from max.genes (top.ranked, top.ranked/2, top.ranked/4, top.ranked/8, top.ranked/16/ top.ranked/32, top.ranked/64). Pearson correlation is calculated on 1000 random selections of genes for each threshold. it is returned the frequence of having pearson >= 0.5.  
#' @examples
#' \dontrun{
#'  library(rCASC)
#'  integrationPsblk(group="docker", 
#'         scratch.folder="/scratch", 
#'         fileX="/data/reanalysis_on_AIsc/comparing_CRC0327/NT_CTX/CRC0327_NT_2_clx/VandE/VandE_bulkRow.csv",
#'         fileY="/data/reanalysis_on_AIsc/comparing_CRC0327/NT_CTX/CRC0327_cetux_2_clx/VandE/VandE_bulkRow.csv", 
#'         separatorX=",",
#'         separatorY=",",
#'         max.genes=320,
#'         outputFolder="/data/reanalysis_on_AIsc/comparing_CRC0327/NT_CTX"
#'  )
#'}
#' @export
integrationPsblk <- function(group=c("sudo","docker"), scratch.folder, fileX, fileY, separatorX, separatorY, max.genes=500, outputFolder){
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


fileX=paste("/scratch/",fileX.file,sep="")
fileY=paste("/scratch/",fileY.file,sep="")

  #executing the docker job

  params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ", data.folder, ":/data -d docker.io/repbioinfo/combinetoprnk Rscript /home/combineExperiments_toprnk.R ", max.genes, " ", fileX, " ", fileY," ", separatorX," ", separatorY, sep="")

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
    system(paste("cp -r ",scrat_tmp.folder,"/XYpb ",data.folder,"/",sep=""))
  
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}
 

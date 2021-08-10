#' @title dePblkae
#' @description This function execute dePblkae analysis which detects differentiale expressed genes among clusters using a pseudobulk generated with autoencoders: autoencoder4pseudoBulk function. It is advisible to  do the analysis on a dataset made of at least 20 permutation
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param data.folder, folder where the total.csv file is located
#' @param scratch.folder, where calculation is performed
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return a file called DE.txt, which contains for each cluster the ID of the genes detected as DE,  a file called log2_reformatedCPM.psblkAE.csv, which contains the log2 CPM of the psblkAE file and a file log2_reformatedCPM.psblkAE_mean-centered.csv, where log2 CPM are mena centered.
#' @examples
#' \dontrun{
#'  library(rCASC)
#'  dePblkae(group="docker", 
#'         data.folder="/somewhere/in/your/PC", 
#'         scratch.folder="/scratch"
#'  )
#'}
#' @export
dePblkae <- function(group=c("sudo","docker"), data.folder, scratch.folder){
 
  outputFolder=data.folder
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

  system(paste("cp ", data.folder, "/total.csv ", scrat_tmp.folder,"/total.csv",sep=""))

  params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -d docker.io/repbioinfo/desc.2021.01 Rscript /home/extractDE.R", sep="")
  resultRun <- runDocker(group=group, params=params)

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  system("rm -fR dockerID")

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

  

  #Copy result folder
  cat("Copying Result Folder")
 
  system(paste("cp -r ",scrat_tmp.folder,"/DEpblkAE ",data.folder,"/",sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}
 

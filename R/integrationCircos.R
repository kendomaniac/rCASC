#' @title integrationCircos
#' @description This function execute toprnk analysis which search for correspondence between clusters of two different experiments using clusters-pseudobulks, zscored on rows, and the cluster specific genes from comet analysis. Thus, the function clustersBulk and cometsc have to be run in the two datasets before their comparison.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param gsea.file, a character string indicating the path of the final_score.csv  generated with gseaXLmHG, file, with file name and extension included. 
#' @param isc.file, a character string indicating the path of the final_score.csv  generated with seuratIntegrationPermutation, file, with file name and extension included. included. 
#' @param XYpb.file, a character string indicating the path of the final_score.csv  generated with toprnk, file, with file name and extension included. included. 
#' @param outputFolder, where results are placed
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return A picture called integrated_score.png and a file called integrated_score.csv and all the final_scores.csv used to produce the integrated results. 
#' @examples
#' \dontrun{
#'  library(rCASC)
#'  integrationCircos(group="docker", 
#'         scratch.folder="/scratch", 
#'         gsea.file="/data/reanalysis_on_AIsc/comparing_CRC0327/NT_CTX/CRC0327_NT_2_clx/VandE/VandE_bulkRow.csv",
#'         isc.file="/data/reanalysis_on_AIsc/comparing_CRC0327/NT_CTX/CRC0327_cetux_2_clx/VandE/VandE_bulkRow.csv",
#'         XYpb.file="/data/reanalysis_on_AIsc/comparing_CRC0327/NT_CTX/CRC0327_cetux_2_clx/VandE/VandE_bulkRow.csv", 
#'         outputFolder="/data/reanalysis_on_AIsc/comparing_CRC0327/NT_CTX"
#'  )
#'}
#' @export
integrationCircos <- function(group=c("sudo","docker"), scratch.folder, gsea.file, isc.file, XYpb.file, outputFolder){
 
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

system(paste("cp ",gsea.file," ",scrat_tmp.folder,"/gsea_final_score.csv",sep=""))
system(paste("cp ",isc.file," ",scrat_tmp.folder,"/isc_final_score.csv",sep=""))
system(paste("cp ",XYpb.file," ",scrat_tmp.folder,"/XYpb_final_score.csv",sep=""))

params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -d docker.io/repbioinfo/xlmhg.2021.01 Rscript /fs_circors.R", sep="")

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
 
    system(paste("cp -r ",scrat_tmp.folder,"/integrated_score ",data.folder,"/",sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}
 

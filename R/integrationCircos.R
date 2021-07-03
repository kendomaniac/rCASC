#' @title integrationCircos
#' @description This function executes integrationCircos, which plots the integrated results of the analysis performed with gseaXLmHG, seuratIntegrationPermutation, integrationPsblk, integrationPblkae and BC
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param gsea.file, a character string indicating the path of the final_score.csv  generated with gseaXLmHG, file, with file name and extension included. 
#' @param isc.file, a character string indicating the path of the final_score.csv  generated with seuratIntegrationPermutation, file, with file name and extension included. included. 
#' @param XYpb.file, a character string indicating the path of the final_score.csv  generated with integrationPsblk, file, with file name and extension included. 
#' @param pblkae.file, a character string indicating the path of the final_score.csv  generated with integrationPblkae, file, with file name and extension included. 
#' @param bcsc.file, a character string indicating the path of the FINAL_score.csv  generated with BCscWrapper, file, with file name and extension included. 
#' @param Xcls.groups, a vector of strings describing the groups of more similar clusters.  The optimal order of the clusters can be deduced by the output of integrationPblkae function with the option type="intra". Format: c("1cl6", "1cl2-1cl3", "1cl1-1cl4-1cl5").
#' @param Ycls.groups, a vector of strings describing the order of the clusters.  The optimal order of the clusters can be deduced by the output of integrationPblkae function with the option type="intra".  Format: c("2cl1-2cl3", "2cl2-2cl4", "2cl5-2cl6")
#' @param outputFolder, where results are placed
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return A picture called integrated_score.png and a file called integrated_score.csv and all the final_scores.csv used to produce the integrated results. The colour ramp for Xgroups is yellow-magenta as instead for the Ygroups is green-blue. Clusters sharing the same colors are those characterized by belonging to the same subgroup detected using integrationPblkae with type="intra". Picture and data used for the integration are lccated in integrated_score folder. Inter cluster edges color code: 1 red; 1-0.7 green; 0.7-0.5 blue; 0.5-0.3 violet; 0.3-0.2 grey; 0.2-0 gold. 
#' @examples
#' \dontrun{
#' library(rCASC)
#' integrationCircos(group="docker", 
#'                   scratch.folder="/scratch", 
#'                   gsea.file="/data/reanalysis_on_AIsc/comparing_CRC0327/NT1_NT2/GSEA/final_score.csv",
#'                   isc.file="/data/reanalysis_on_AIsc/comparing_CRC0327/NT1_NT2/ISC/final_score.csv",
#'                   XYpb.file="/data/reanalysis_on_AIsc/comparing_CRC0327/NT1_NT2/XYpb/XYpb_final_score.csv",
#'                   pblkae.file="/data/reanalysis_on_AIsc/comparing_CRC0327/NT1_NT2/pblkAE/final_score.csv",
#'                   bcsc.file="/data/reanalysis_on_AIsc/comparing_CRC0327/NT1_NT2/BCsc/FINAL_score.csv",
#'                   Xcls.order=NULL,
#'                   Ycls.order=NULL, 
#'                   outputFolder="/data/reanalysis_on_AIsc/comparing_CRC0327/NT1_NT2"
#' )
#' 
#' integrationCircos(group="docker", 
#'                   scratch.folder="/scratch", 
#'                   gsea.file="/data/reanalysis_on_AIsc/comparing_CRC0327/NT1_NT2/GSEA/final_score.csv",
#'                   isc.file="/data/reanalysis_on_AIsc/comparing_CRC0327/NT1_NT2/ISC/final_score.csv",
#'                   XYpb.file="/data/reanalysis_on_AIsc/comparing_CRC0327/NT1_NT2/XYpb/XYpb_final_score.csv",
#'                   pblkae.file="/data/reanalysis_on_AIsc/comparing_CRC0327/NT1_NT2/pblkAE/final_score.csv",
#'                   bcsc.file="/data/reanalysis_on_AIsc/comparing_CRC0327/NT1_NT2/BCsc/FINAL_score.csv",
#'                   Xcls.groups=c("1cl6", "1cl2-1cl3", "1cl1-1cl4-1cl5"),
#'                   Ycls.groups=c("2cl1-2cl3", "2cl2-2cl4", "2cl5-2cl6"), 
#'                   outputFolder="/data/reanalysis_on_AIsc/comparing_CRC0327/NT1_NT2"
#' )
#' 
#'}
#' @export
integrationCircos <- function(group=c("sudo","docker"), scratch.folder, gsea.file=NULL, isc.file=NULL, XYpb.file=NULL, pblkae.file=NULL, bcsc.file=NULL, Xcls.groups=NULL, Ycls.groups=NULL, outputFolder){
 
  if(!is.null(Xcls.groups)){
    Xcls.groups <- paste(Xcls.groups, collapse="_")
  }
  if(!is.null(Ycls.groups)){
    Ycls.groups <- paste(Ycls.groups, collapse="_")
  }
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
  if(!is.null(gsea.file)){
    if(file.exists(gsea.file)){
      system(paste("cp ",gsea.file," ",scrat_tmp.folder,"/gsea_final_score.csv",sep=""))
    }
  }
  if(!is.null(isc.file)){
    if(file.exists(isc.file)){
      system(paste("cp ",isc.file," ",scrat_tmp.folder,"/isc_final_score.csv",sep=""))
    }
  }
  if(!is.null(XYpb.file)){
    if(file.exists(XYpb.file)){
      system(paste("cp ",XYpb.file," ",scrat_tmp.folder,"/XYpb_final_score.csv",sep=""))
    }
  }
  if(!is.null(pblkae.file)){
    if(file.exists(pblkae.file)){
      system(paste("cp ",pblkae.file," ",scrat_tmp.folder,"/pblkae_final_score.csv",sep=""))
    }
  }
  if(!is.null(bcsc.file)){
    if(file.exists(bcsc.file)){
      system(paste("cp ",bcsc.file," ",scrat_tmp.folder,"/bcsc_final_score.csv",sep=""))
    }
  }
  
if(is.null(Xcls.groups) && is.null(Ycls.groups)){
  params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -d docker.io/repbioinfo/xlmhg.2021.01 Rscript /home/fs_circors.R", sep="")
}else{
  params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -d docker.io/repbioinfo/xlmhg.2021.01 Rscript /home/fs_circors1.R ", Xcls.groups, " ", Ycls.groups, sep="")
}


resultRun <- runDocker(group=group, params=params)

  #waiting for the end of the container work
  if(resultRun==0){
    #system(paste("cp  ", scrat_tmp.folder, "/*top* ", data.folder, sep=""))
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
 

#' @title integrationPblkae
#' @description This function execute integrationPblkae analysis which search for correspondence between clusters of two different experiments using clusters-pseudobulks generated using sparsely connected autoencoders. Thus, the function autoencoder4pseudoBulk has to be run in the two datasets before their comparison.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param fileX, a character string indicating the path to the total.csv.for the 1st dataset to be integrated. Total.csv is generated with autoencoder4pseudoBulk. File, with file name and extension included. 
#' @param fileY, a character string indicating the path to the total.csv.for the 2nd dataset to be integrated. Total.csv is generated with autoencoder4pseudoBulk. File, with file name and extension included. 
#' @param type, two values inter, intra. Inter refers to comparison among clusters of two independent experiments. Intra comparisons among clusters of the same experiment
#' @param stats, two values anovalike, pairwise. Anovalike refers to comparison among clusters with respect to a pseudo reference sample. pairwise comparisons among all clusters by pairwise comparison
#' @param outputFolder, where results are placed
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return a folder called psblkAE, which contains file called final_score.csv and all the intermediate files used to produce the integrated results. 
#' @examples
#' \dontrun{
#'  #inter
#'  library(rCASC)
#'  integrationPblkae(group="docker", 
#'         scratch.folder="/scratch", 
#'         fileX="/data/clusters_association_paper/setA1_set1/setA1/VandE/Results/setA1/permutation/total.csv",
#'         fileY="/data/clusters_association_paper/setA1_set1/set1/VandE/Results/set1/permutation/total.csv",
#'         outputFolder="/data/clusters_association_paper/setA1_set1",
#'         type="inter"
#'         stats="anovalike"
#'  )
#'  
#'  #intra
#'  library(rCASC)
#'  integrationPblkae(group="docker", 
#'         scratch.folder="/scratch", 
#'         fileX="/data/clusters_association_paper/setA1_set1/setA1/VandE/Results/setA1/permutation/total.csv",
#'         outputFolder="/data/clusters_association_paper/setA1_set1/setA1",
#'         type="intra"
#'         stats="pairwise"
#'  )
#'}
#' @export
integrationPblkae <- function(group=c("sudo","docker"), scratch.folder, fileX=NULL, fileY=NULL, outputFolder, type=c("inter", "intra"), stats=c("anovalike", "pairwise")){
 
  if(type=="intra"){
    fileY=fileX
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

system(paste("cp ", fileX," ",scrat_tmp.folder,"/X_total.csv",sep=""))
system(paste("cp ", fileY," ",scrat_tmp.folder,"/Y_total.csv",sep=""))

if(stats=="anovalike"){
	params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -d docker.io/repbioinfo/desc.2021.01 Rscript /home/build.anova.like.dataset.R", sep="")
	resultRun <- runDocker(group=group, params=params)
}else if(stats=="pairwise"){
	params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -d docker.io/repbioinfo/desc.2021.01 Rscript /home/build.pairwise.dataset.R", sep="")
	resultRun <- runDocker(group=group, params=params)
}

#ending the first part
cat("\npreprcessing of clusters of fileX and fileY done\n")

#saving log and removing docker container
container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/", substr(container.id,1,12),".log", sep=""))
system(paste("docker rm ", container.id, sep=""))
system("rm -fR dockerID")

if(stats=="anovalike"){
	params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -d docker.io/repbioinfo/desc.2021.01 Rscript /home/debulkAE.R /scratch/X_counts_reformat.txt X ", type, sep="")
	resultRun <- runDocker(group=group, params=params)

	cat("\nanovaLike for fileX done\n")

	#saving log and removing docker container
	container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
	system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/", substr(container.id,1,12),".log", sep=""))
	system(paste("docker rm ", container.id, sep=""))
	system("rm -fR dockerID")

	if(type=="inter"){
  	  	params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -d docker.io/repbioinfo/desc.2021.01 Rscript /home/debulkAE.R /scratch/Y_counts_reformat.txt Y ", type, sep="")
  		resultRun <- runDocker(group=group, params=params)
  
  		cat("\nanovaLike for fileY done\n")
  
 	 	#saving log and removing docker container
  		container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  		system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/", substr(container.id,1,12),".log", sep=""))
  		system(paste("docker rm ", container.id, sep=""))
  		system("rm -fR dockerID")
	}


	params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -d docker.io/repbioinfo/desc.2021.01 Rscript /home/post_debulkAE.R ", type, sep="")
	resultRun <- runDocker(group=group, params=params)

	cat("\nintegration done\n")

	container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
	system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/", substr(container.id,1,12),".log", sep=""))
	system(paste("docker rm ", container.id, sep=""))
	system("rm -fR dockerID")
	system("rm  -fR tempFolderID")
}else if(stats=="pairwise"){
	params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -d docker.io/repbioinfo/desc.2021.01 Rscript debulkAE_pairwise ", type, sep="")
	resultRun <- runDocker(group=group, params=params)

	cat("\nanovaLike for fileX done\n")

	#saving log and removing docker container
	container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
	system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/", substr(container.id,1,12),".log", sep=""))
	system(paste("docker rm ", container.id, sep=""))
	system("rm -fR dockerID")
	
	
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

  

  #Copy result folder
  cat("Copying Result Folder")
 
  system(paste("cp -r ",scrat_tmp.folder,"/pblkAE ",data.folder,"/",sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}
 

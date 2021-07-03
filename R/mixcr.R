#' @title mixcr
#' @description This function creates mixcrFiles from bulk TCRseq
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param fastqPath, Path of fastq folder
#' @param resFolderCustom, optional parameter. Default will store the results in fastqPath otherwise will store the results in resFolderCustom path. 
#' @author Luca Alessandrì
#'
#'
#' @return a formatted mixcr file
#' @examples
#' \dontrun{
#' library(rCASC)
#' dir.create("scratch")
#' scratch.folder=paste(getwd(),"scratch",sep="/")
#' fastqPath=paste(getwd(),"fastq",sep="/")
#' resFolder=paste(getwd(),"resFolder",sep="/")
#' dir.create(resFolder)
#' mixcr(group="docker",scratch.folder=scratch.folder,fastqPath=fastqPath,resFolderCustom=resFolder)
#' }
#'
#'
#' @export

mixcr <- function(group=c("sudo","docker"),scratch.folder,fastqPath,resFolderCustom="NULL"){


genomeFolder=fastqPath
      dockerImage="repbioinfo/mixcr:1"
  
    setwd(genomeFolder)

#storing the position of the home folder
  home <- getwd()


  #running time 1
  ptm <- proc.time()

  #setting the data.folder as working folder


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



  #check  if scratch folder exist
  if (!file.exists(scratch.folder)){
    cat(paste("\nIt seems that the ",scratch.folder, " folder does not exist\n"))
    system("echo 3 > ExitStatusFile 2>&1")
      setwd(home)
    return(3)
  }
  tmp.folder <- gsub(":","-",gsub(" ","-",date()))
  scrat_tmp.folder=file.path(scratch.folder, tmp.folder)
  writeLines(scrat_tmp.folder,paste(genomeFolder,"/tempFolderID", sep=""))
  cat("\nCreating a folder in scratch folder\n")
  scrat_tmp.folder=file.path(scrat_tmp.folder)
  dir.create(scrat_tmp.folder)

  #cp fastq folder in the scrat_tmp.folder
  #executing the docker job

  params <- paste("--cidfile ",genomeFolder,"/dockerID  -v ", scrat_tmp.folder, ":/scratch -d ",dockerImage, " /home/mixcr2.sh ",basename(genomeFolder),sep="")

system(paste("cp -r ",genomeFolder,"/* ",scrat_tmp.folder,sep=""))
  #Run docker
  resultRun <- runDocker(group=group, params=params)
  #waiting for the end of the container work
  if(resultRun==0){
  if(resFolderCustom=="NULL"){
  res=paste(genomeFolder,"Results",sep="/")}
  else{res=resFolderCustom}
    dir.create(res)
    system(paste("cp -R ", scrat_tmp.folder, "/*.txt ",res, sep=""))
        system(paste("cp -R ", scrat_tmp.folder, "/*.fasta ",res, sep=""))
                system(paste("cp -R ", scrat_tmp.folder, "/trim ",res, sep=""))


  }

 #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(genomeFolder)
  dir <- dir[grep("run.info",dir)]
  if(length(dir)>0){
    con <- file("run.info", "r")
    tmp.run <- readLines(con)
    close(con)
    tmp.run[length(tmp.run)+1] <- paste("cellranger user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("cellranger system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("cellranger elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
     tmp.run[1] <- paste("cellranger run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("cellranger system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("cellranger elapsed run time mins ",ptm[3]/60, sep="")

    #writeLines(tmp.run,"run.info")
  }

  #saving log and removing docker container
  container.id <- readLines(paste(genomeFolder,"/dockerID", sep=""), warn = FALSE)
  #system(paste("docker logs ", substr(container.id,1,12), " &> ",genomeFolder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  #removing temporary folder
#  cat("\n\nRemoving the temporary file ....\n")
#  system(paste("rm -fR ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system("rm -fR ExitStatusFile")

  #system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",genomeFolder, sep=""))
   system(paste("rm  -fR ",res,"/tempFolderID",sep=""))

   system(paste("rm -fR ",res,"/ExitStatusFile",sep=""))
   system(paste("chmod -R 777 ",res,sep=""))


 
 
 
 
 
 
 setwd(home)

} 

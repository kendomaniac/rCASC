#' @title Marker genes discovery with COMETSC
#' @description This function executes a ubuntu docker for cometsc (https://github.com/MSingerLab/COMETSC)
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param file, a character string indicating the path of the file, with file name and extension included
#' @param threads, integer refering to the max number of process run in parallel default 1 max the number of clusters under analysis, i.e. nCluster
#' @param X, from 0 to 1 argument for XL-mHG default 0.15, for more info see cometsc help.
#' @param K, the number of gene combinations to be considered., possible values 2, 3, 4, default 2. WARNING increasing the number of combinations makes the matrices very big
#' @param counts, If set to True it will graph the log(expression+1). To be used if unlogged data are provided
#' @param skipvis, Set to True to skip visualizations
#' @param nCluster, number of interested cluster used for analysis
#' @param scratch.folder, temporary folder where calculation is made
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @return folders with prefix output. More info in output at https://hgmd.readthedocs.io/en/latest/Output.html
#' @author Raffaele Calogero,raffaele.calogero [at] unito [dot] it, University of Torino
#' 
#' @examples
#' \dontrun{
#'     #running cometsc
#'     cometsc(group="docker", file="/Users/raffaelecalogero/Desktop/AXLN1/data/topx_veanno.csv", 
#'            scratch.folder="/Users/raffaelecalogero/Desktop",
#'            threads=1, counts="True", skipvis="False", nCluster=8, separator=",") 
#' }
#'
#' @export
cometsc <- function(group=c("sudo","docker"), file, scratch.folder, threads=1,  X=0.15, K=2, counts=c("True", "False"), skipvis=c("True", "False"), nCluster, separator){

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
  system(paste("cp -r ",data.folder,"/Results/", matrixName,"/",nCluster, "/* ",scrat_tmp.folder,sep=""))
  system(paste("cp -r ", file," ",scrat_tmp.folder,sep=""))
  
  #executing the docker job
  params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ", data.folder, ":/data -d docker.io/repbioinfo/cometsc.2020.01 sh /bin/cometsc.sh ", matrixName, " ", threads, " ", X, " ", K, " ", counts, " ", skipvis, " ", nCluster," ", separator, sep="")
  resultRun <- runDocker(group=group, params=params)
  
  #waiting for the end of the container work
  if(resultRun==0){
    system(paste("cp -r ", scrat_tmp.folder,"/output* ",data.folder,"/Results/", sep=""))
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
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="docker4seq"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}

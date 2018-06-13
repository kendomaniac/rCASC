#' @title A function to handle sigle cell Lorenz Quality filter for Single-cells
#' @description This function executes a docker that produces ....
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param data.folder, a character string indicating the folder where input data are located and where output will be written
#' @param matrixName, counts table name. Matrix data file must be in data.folder. The file MUST contain RAW counts, without any modification, such as log transformation, normalizatio etc. 
#' @param p_value, threshold to be used for the filtering
#' @param format, counts file extension, "txt", "csv"
#' @param separator, separator used in count file, e.g. '\\t', ','
#' 
#' @author Name Family name, myemail [at] somewhere [dot] org, Affiliation
#'
#' @return output will be in the same format and with the same separator of input.
#' 
#' @examples
#' \dontrun{
#'         system("wget http://130.192.119.59/public/testSCumi_mm10.csv.zip")
#'      library(casc)
#'      system("unzip testSCumi_mm10.csv.zip")
#'      #filtering low quality cells
#'      lorenzFilter(group="docker",scratch.folder="/data/scratch/",
#'                   file=paste(getwd(),"testSCumi_mm10.csv",sep="/"),
#'                   p_value=0.05,separator=',')
#' 
#' @export
lorenzFilter <- function(group=c("sudo","docker"), scratch.folder, file, p_value, separator){

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

   write.csv(read.table(paste(data.folder,"/",matrixName,".",format,sep=""),sep=separator,header=TRUE,row.names=1),paste(scrat_tmp.folder,"/set1.csv",sep=""))
if(separator=="\t"){
separator="tab"
}

  #executing the docker job
      params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ", data.folder, ":/data -d docker.io/rcaloger/lorenz Rscript /home/main.R ",matrixName," ",p_value," ",format," ",separator, sep="")
 
resultRun <- runDocker(group=group, params=params)
  
  #waiting for the end of the container work
  if(resultRun==0){
  #  system(paste("cp ", scrat_tmp.folder, "/* ", data.folder, sep=""))
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
  system(paste("cp -r ",scrat_tmp.folder,"/* ",data.folder,"/Results/",sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  #system(paste("cp ",paste(path.package(package="docker4seq"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}

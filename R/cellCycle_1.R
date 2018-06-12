#' @title Cell Cycle 
#' @description This function executes a ubuntu docker that associates to each cell a cell cycle stage
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param data.folder, a character string indicating the folder where input data are located and where output will be written
#' @param matrixName, counts table name. Matrix data file must be in data.folder. The file MUST contain RAW counts, without any modification, such as log transformation, normalizatio etc. 
#' @param format, matrix count format, "csv", "txt"#' @param B, second Cluster that has to be merged
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param geneNameControl, 0 if the matrix has gene name without ENSEMBL code.1 if the gene names is formatted like this : ENSMUSG00000000001:Gnai3. If the gene names is only ensamble name you have to run SCannoByGtf before start using this script.
#' @param window, this parameters let you see the cell cycle algorithm plot in block.  
#' @param seed, important parameter for reproduce the same result with the same input
#' @author Luca Alessandri , alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return will change all the files generated from permAnalysis algorithm in a new folder matrixName_Cluster_merged/
#' @examples
#'\dontrun{
#' cellCycle(group,scratch.folder,data.folder,mainMatrix,format,separator,window,seed)
#'}
#' @export

cellCycle1 <- function(group=c("sudo","docker"), scratch.folder, file,separator,geneNameControl=0,window=1,seed=111){

  data.folder=dirname(file)
positions=length(strsplit(basename(file),"\\.")[[1]])
matrixNameC=strsplit(basename(a),"\\.")[[1]]
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
if(separator=="\t"){
separator2="tab"
}else{separator2=separator}

#gene Name Control
if(geneNameControl==0){
system(paste("cp ",data.folder,"/",matrixName,".",format," ",data.folder,"/",matrixName,"_old.",format,sep=""))
mainMatrix=read.table(paste(data.folder,"/",matrixName,".",format,sep=""),header=TRUE,row.names=1,sep=separator)
rownames(mainMatrix)=paste(seq(1,nrow(mainMatrix)),":",rownames(mainMatrix),sep="")
write.table(mainMatrix,paste(data.folder,"/",matrixName,".",format,sep=""),sep=separator,col.names=NA)
}



system(paste("cp ",data.folder,"/",matrixName,".",format," ",scrat_tmp.folder,sep=""))


  
  #executing the docker job
params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ", data.folder, ":/data -d docker.io/rcaloger/cellcycle_1 Rscript /home/main.R ",matrixName," ",format," ",separator2," ",window," ",seed,sep="")
resultRun <- runDocker(group=group, params=params)
  
  #waiting for the end of the container work
  if(resultRun==0){
    #system(paste("cp ", scrat_tmp.folder, "/* ", data.folder, sep=""))
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
    dir.create(paste(data.folder,"/Results",sep=""))
        dir.create(paste(data.folder,"/Results/",matrixName,sep=""))

  system(paste("cp ",scrat_tmp.folder,"/cellCycleRange.pdf ",data.folder,"/Results/",matrixName,sep=""))

  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
#  system(paste("cp ",paste(path.package(package="docker4seq"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}

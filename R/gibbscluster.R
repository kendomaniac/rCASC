#' @title gibbsR
#' @description This function creates performs motif search in peptides. It was developed to use gibbscluster on TCR CDR3 peptide sequences
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, path to peptide file, 
#' @param resFolderCustom, optional parameter. Default will store the results in fastqPath otherwise will store the results in resFolderCustom path. 
#' @param jobName, identify of the sample
#' @param nCluster, number of clusters to be generated
#' @param motifLength, length of the motif to be searched
#' @param maxDelLength, max length of deletions in the peptide search, expressed in number of aa
#' @param maxInsLength, max length of insertion in the peptide search, expressed in number of aa
#' @param numbOfSeed, number of initial configurations to be used to detect the highest KLD sum
#' @param penalityFactorIntCluster, default 0.8,
#' @param backGroundAminoFreq, two value 1 defined by uniprot and 2 defined by the dataset under analysis
#' @param seqWeightType, sequence weighting type: default 1, faster but less precise, 2 slower but more precise, 3 none
#' @author Luca Alessandrì
#'
#'
#' @return an indexed genome compliant with 10XGenomics cellranger
#' @examples
#' \dontrun{
#' library(rCASC)
#' dir.create("scratch")
#' scratch.folder=paste(getwd(),"scratch",sep="/")
#' fastqPath=paste(getwd(),"fastq",sep="/")
#' resFolder=paste(getwd(),"resFolder",sep="/")
#' dir.create(resFolder)
#' 
#' gibbsR(group="docker",scratch.folder=scratch.folder,file=file,resFolderCustom=results,jobName="test",nCluster=8,motifLength=4,maxDelLength=4,maxInsLength=4,numbOfSeed=4,penalityFactorIntCluster=0.8,backGroundAminoFreq=2,seqWeightType=1) 
#' 
#' newFolder=paste(results,"res",list.files(paste(results,"res",sep="/")),"cores",sep="/")
#' }
#'
#'
#' @export



gibbsR <- function(group=c("sudo","docker"),scratch.folder,file,resFolderCustom="NULL",jobName,nCluster,motifLength,maxDelLength,maxInsLength,numbOfSeed,penalityFactorIntCluster=0.8,backGroundAminoFreq=2,seqWeightType=1){
  
  name=basename(file)
  genomeFolder=dirname(file)
  dockerImage="repbioinfo/mixcr:v1"
  
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
  
  params <- paste("--cidfile ",genomeFolder,"/dockerID  -v ", scrat_tmp.folder, ":/scratch -d ",dockerImage, " /home/gscript.sh ",name," ",jobName," ",nCluster," ",motifLength," ",maxDelLength," ",maxInsLength," ",numbOfSeed," ",penalityFactorIntCluster," ",backGroundAminoFreq," ",seqWeightType,sep="")
  system(paste("cp ",file," ",scrat_tmp.folder,sep=""))
  #Run docker
  resultRun <- runDocker(group=group, params=params)
  
  #waiting for the end of the container work
  
  if(resFolderCustom=="NULL"){
    res=paste(genomeFolder,"Results",sep="/")}
  else{res=resFolderCustom}
  dir.create(res)
  system(paste("cp -r ", scrat_tmp.folder, "/res/"," ",res, sep=""))
  
  
  
  #saving log and removing docker container
  container.id <- readLines(paste(genomeFolder,"/dockerID", sep=""), warn = FALSE)
  #system(paste("docker logs ", substr(container.id,1,12), " &> ",genomeFolder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  #removing temporary folder
  #  cat("\n\nRemoving the temporary file ....\n")
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system("rm -fR ExitStatusFile")
  
  #system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",genomeFolder, sep=""))
  system(paste("rm  -fR ",res,"/tempFolderID",sep=""))
  
  system(paste("rm -fR ",res,"/ExitStatusFile",sep=""))
  system(paste("chmod -R 777 ",res,sep=""))
  #system(paste("rm -fR ",scrat_tmp.folder))
  setwd(home)
  
} 
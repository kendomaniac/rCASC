#' @title Spatial transcriptomics pipeline
#' @description Create count matrix from spatial transcriptomics fasta
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param data.folder, a character string indicating the path of the result folder.
#' @param genome.folder, a character string indicating the path of the genome folder. 
#' @param fastqPathFolder, a character string indicating the path of fastq folder
#' @param ID, a character string indicating the name of the project
#' @param imgNameAndPath, path and name of tiff image file required for analysis. 
#' @param slide, identificative number from dataset download
#' @param area, identificative value from dataset download 

#' @author Luca Alessandri , alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return count matrix from spatial transcriptomics
#' @examples
#'\dontrun{
#' Dataset="curl -O http://s3-us-west-2.amazonaws.com/10x.files/samples/spatial-exp/1.0.0/V1_Mouse_Kidney/V1_Mouse_Kidney_fastqs.tar"
#' DatasetImage="curl -O http://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Mouse_Kidney/V1_Mouse_Kidney_image.tif"
#' referenceGenomeHG38="curl -O http://cf.10xgenomics.com/supp/spatial-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz"
#' referenceGenomeMM10="curl -O http://cf.10xgenomics.com/supp/spatial-exp/refdata-cellranger-mm10-3.0.0.tar.gz"
#' stpipeline(group="docker", scratch.folder="/run/media/user/Maxtor4/scratch", data.folder="/run/media/user/Maxtor4/prova2", genome.folder="/home/user/spatial/refdata-cellranger-mm10-3.0.0", fastqPathFolder="/home/user/spatial/V1_Mouse_Kidney_fastqs", ID="hey",imgNameAndPath="/home/user/spatial/V1_Mouse_Kidney_image.tif",slide="V19L29-096",area="B1")
#'}
#' @export
stpipeline <- function(group=c("sudo","docker"), scratch.folder,data.folder,genome.folder,fastqPathFolder,ID,imgNameAndPath,slide=NULL,area=NULL){


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

  #preprocess matrix and copying files



dir.create(paste(scrat_tmp.folder))
  dir.create(file.path(paste(scrat_tmp.folder,"/fastq",sep="")))
system(paste("cp -r ",fastqPathFolder,"/* ",scrat_tmp.folder,"/fastq",sep=""))
system(paste("cp -r ",imgNameAndPath," ",scrat_tmp.folder,"/",sep=""))
 imgname=basename(imgNameAndPath)

  #executing the docker job
    #params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ",genome.folder,":/genome -d docker.io/rcaloger/stpipeline python3 /data2/st_pipeline-master/scripts/st_pipeline_run.py --expName ",nameExp," --ids /scratch/",ids," --ref-map /genome/ --log-file log_file.txt --output-folder /data --ref-annotation /scratch/",gtf," /scratch/",f1," /scratch/",f2,sep="")
if(is.null(slide) & is.null(area)){
  params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ",genome.folder,":/genome -d repbioinfo/stpipelineofficial /home/spatial.sh ",ID," ",basename(imgNameAndPath),sep="")
}else{
  params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ",genome.folder,":/genome -d repbioinfo/stpipelineofficial /home/spatial2.sh ",ID," ",basename(imgNameAndPath)," ",slide," ",area,sep="")

}
resultRun <- runDocker(group=group, params=params)

 
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
  system(paste("cp -r ",scrat_tmp.folder,"/* ",data.folder,sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
} 
 

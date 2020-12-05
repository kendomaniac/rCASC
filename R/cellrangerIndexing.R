#' @title Cellranger indexing
#' @description This function creates the indexing for 10Xgenomics
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param gtf.url, a character string indicating the URL from ENSEMBL ftp for the GTF for genome of interest
#' @param fasta.url, a character string indicating the URL from ENSEMBL ftp for the unmasked genome sequence of interest
#' @param genomeFolder,  path for the genome folder
#' @param bio.type, ENSEMBL biotype to filter the GTF
#' @param nThreads, number of cores for parallelization
#' @param version,  cellranger version: 2, 3 or 5. 
#' @author Luca Alessandrì
#'
#'
#' @return an indexed genome compliant with 10XGenomics cellranger
#' @examples
#' \dontrun{
#' library(rCASC)
#' setwd("/data/genomes/hg38refcellranger")
#'
#' cellrangerIndexing(group="docker", scratch.folder="/data/scratch", 
#'             gtf.url="ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz",
#'             fasta.url="ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
#'             genomeFolder = getwd(), bio.type="protein_coding", nThreads = 8)
#' }
#'
#'
#' @export

cellrangerIndexing <- function(group=c("sudo","docker"),scratch.folder,genomeFolder,gtf.url,fasta.url,bio.type=c("protein_coding","unitary_pseudogene",
                                                           "unprocessed_pseudogene","processed_pseudogene",
                                                           "transcribed_unprocessed_pseudogene","processed_transcript",
                                                           "antisense","transcribed_unitary_pseudogene",
                                                           "polymorphic_pseudogene","lincRNA",
                                                           "sense_intronic","transcribed_processed_pseudogene",
                                                           "sense_overlapping","IG_V_pseudogene",
                                                           "pseudogene","TR_V_gene",
                                                           "3prime_overlapping_ncRNA","IG_V_gene",
                                                           "bidirectional_promoter_lncRNA","snRNA",
                                                           "miRNA","misc_RNA",
                                                           "snoRNA","rRNA",
                                                           "IG_C_gene","IG_J_gene",
                                                           "TR_J_gene","TR_C_gene",
                                                           "TR_V_pseudogene","TR_J_pseudogene",
                                                           "IG_D_gene","ribozyme",
                                                           "IG_C_pseudogene","TR_D_gene",
                                                           "TEC","IG_J_pseudogene",
                                                           "scRNA","scaRNA",
                                                           "vaultRNA","sRNA",
                                                           "macro_lncRNA","non_coding","IG_pseudogene"),nThreads, version="5"){

  id="results_cellranger"
  #docker image

 
    if(version == "2"){
      dockerImage="docker.io/repbioinfo/cellranger"
    } else if(version == "3"){
      dockerImage="docker.io/repbioinfo/cellranger.2018.03"
    } else if(version == "5"){
      dockerImage="docker.io/repbioinfo/cellranger.2020.05"
    }
    
  

#storing the position of the home folder
  home <- getwd()


  #running time 1
  ptm <- proc.time()

  #setting the data.folder as working folder
  if (!file.exists(genomeFolder)){
    cat(paste("\nIt seems that the ",genomeFolder, " folder does not exist\n"))
    return(2)
  }
  setwd(genomeFolder)

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

  params <- paste("--cidfile ",genomeFolder,"/dockerID  -v ", scrat_tmp.folder, ":/data -d ",dockerImage, " /home/indexing.sh ",gtf.url," ",fasta.url," ",bio.type," ",nThreads,sep="")

#system(paste("cp -r ",genomeFolder,"/* ",scrat_tmp.folder,sep=""))
  #Run docker
  resultRun <- runDocker(group=group, params=params)
  #waiting for the end of the container work
  if(resultRun==0){
    system(paste("cp -R ", scrat_tmp.folder, "/* ",genomeFolder, sep=""))
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

    writeLines(tmp.run,"run.info")
  }

  #saving log and removing docker container
  container.id <- readLines(paste(genomeFolder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", substr(container.id,1,12), " &> ",genomeFolder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  #removing temporary folder
#  cat("\n\nRemoving the temporary file ....\n")
#  system(paste("rm -fR ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",genomeFolder, sep=""))
 
 

 
 
 
 
 
 
 setwd(home)

} 

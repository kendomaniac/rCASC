#' @title Cellranger indexing
#' @description This function creates the indexing for 10Xgenomics
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param gtf.url, url for gtf download
#' @param fasta.url, url for fasta download
#' @param genomeFolder,  path for the genome folder
#' @param bio.type, biotype to filter the gtf
#' @params nThreads, number of cores for parallelization 
#' @author Luca Alessandrì
#'
#'
#' @return a folder called results_cellranger, more info on the structure of this folder at https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview . In /somewhewre_in_your_computer/results_cellranger/outs/filtered_gene_bc_matrices the cells counts matrices results_cellranger.cvs and results_cellranger.txt are saved for further use.
#'
#' @examples
#' \dontrun{
#' home <- getwd()
#' library(rCASC)
#' setwd("/data/genomes/cellranger_hg19mm10")
#' #getting the human and mouse cellranger index
#' system("wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-hg19-and-mm10-2.1.0.tar.gz")
#' untar("refdata-cellranger-hg19-and-mm10-2.1.0.tar.gz") 
#' setwd(home)
#' # 100 cells 1:1 Mixture of Fresh Frozen Human (HEK293T) and Mouse (NIH3T3) Cells
#' system("wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/hgmm_100/hgmm_100_fastqs.tar")
#' untar("hgmm_100_fastqs.tar")
#' home=paste(home,"/fastqs",sep="")
#' cellrangerCount(group="docker",  transcriptome.folder="/data/genomes/cellranger_hg19mm10/refdata-cellranger-hg19_and_mm10-2.1.0",  fastq.folder=getwd(),  expect.cells=100, nosecondary=TRUE, scratch.folder="/data/scratch", version="2")
#' 
#' sraDownload(group = "docker", sra.name = "SRR7762358", data.folder = getwd(), scratch.folder = "/data/scratch", threads = 8)
#' system("mv ./SRR7762358/SRR7762358.fastq.gz ./SRR7762358/SRR7762358_S1_L001_R1_001.fastq.gz")
#' cellrangerCount(group="docker",  transcriptome.folder="/data/genomes/refdata-cellranger-GRCh38-3.0.0",  fastq.folder=getwd(), sample="SRR7762358",  nosecondary=TRUE, scratch.folder="/data/scratch", version="3")
#'
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
                                                           "macro_lncRNA","non_coding","IG_pseudogene"),nThreads){

  id="results_cellranger"
  #docker image

 
    dockerImage="docker.io/repbioinfo/cellranger.2018.03"
  
  

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

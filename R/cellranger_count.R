#' @title Cellranger count
#' @description This function takes FASTQ files from cellranger mkfastq and performs alignment, filtering, barcode counting, and UMI counting.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param transcriptome.folder,  path to the Cell Ranger compatible transcriptome reference e.g. for a human and mouse mixture sample, use refdata-cellranger-hg19-and-mm10-1.2.0
#' @param fastq.folder,  path of the fastq path folder in fastq folder the  fastq must have the format SAMPLENAME_S1_L001_R1_001.fastq.gz
#' @param sample, fastq name, if fastq name is subject1_S1_L001_R1_001.fastq.gz sample is subject1
#' @param expect.cells,  optional setting the number of recovered cells. Default: 3000 cells.
#' @param force.cells,  optional to force pipeline to use this number of cells, bypassing the cell detection algorithm. Use this if the number of cells estimated by Cell Ranger is not consistent with the barcode rank plot.
#' @param nosecondary,  optional flag, default TRUE, to skip secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Set this if you plan to use cellranger reanalyze or your own custom analysis.
#' @param chemistry,  optional assay configuration. One of: auto for autodetection (default), threeprime for Single Cell 3end, fiveprime for Single Cell 5end, SC3Pv1 for Single Cell 3end v1, SC3Pv2 for Single Cell 3end v2, SC5P-PE for Single Cell 5end paired-end (both R1 and R2 are used for alignment), SC5P-R2 for Single Cell 5end R2-only (where only R2 is used for alignment).
#' @param r1.length,  optional hard-trim the input R1 sequence to this length. Note that the length includes the Barcode and UMI sequences so do not set this below 26 for Single Cell 3end v2 or Single Cell 5end. This and --r2-length are useful for determining the optimal read length for sequencing.
#' @param r2.length,  optional hard-trim the input R2 sequence to this length.
#' @param lanes,  optional, lanes associated with this sample
#' @param localcores,  restricts cellranger to use specified number of cores to execute pipeline stages. By default, cellranger will use all of the cores available on your system.
#' @param localmem,  restricts cellranger to use specified amount of memory, in GB, to execute pipeline stages. By default, cellranger will use 90\% of the memory available on your system. Please note that cellranger requires at least 16 GB of memory to run all pipeline stages.
#' @param version,  cellranger version: 2 or 3. version 3 handle better the counts table generation.
#' @author Greta Romano, romano [dot] greta [at] gmail [dot] com, University of Torino
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

cellrangerCount <- function(group=c("sudo","docker"),  transcriptome.folder,  fastq.folder,  sample=NULL, expect.cells=NULL, force.cells=NULL, nosecondary=TRUE, chemistry=NULL, r1.length=NULL,  r2.length=NULL, lanes=NULL, localcores=NULL, localmem=NULL,  scratch.folder, version="3"){

  id="results_cellranger"
  #docker image
  if(version == "2"){
    dockerImage="docker.io/grromano/cellranger"
  } else if(version == "3"){
    dockerImage="docker.io/repbioinfo/cellranger.2018.03"
  }
  

#storing the position of the home folder
  home <- getwd()


  #running time 1
  ptm <- proc.time()

  #setting the data.folder as working folder
  if (!file.exists(fastq.folder)){
    cat(paste("\nIt seems that the ",fastq.folder, " folder does not exist\n"))
    return(2)
  }
  setwd(fastq.folder)

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
  writeLines(scrat_tmp.folder,paste(fastq.folder,"/tempFolderID", sep=""))
  cat("\nCreating a folder in scratch folder\n")
  scrat_tmp.folder=file.path(scrat_tmp.folder)
  dir.create(scrat_tmp.folder)

  #cp fastq folder in the scrat_tmp.folder
  system(paste("cp ", fastq.folder, "/*.gz ", scrat_tmp.folder, sep=""))

  #executing the docker job
  params <- paste("--cidfile ",fastq.folder,"/dockerID -v ",transcriptome.folder,":/transcr -v ", scrat_tmp.folder, ":/data -d ",dockerImage, " /bin/cellranger count  --id=",id," --transcriptome=/transcr --fastqs=/data", sep="")

  if(!is.null(sample)){

    params<-paste(params," --sample=",sample, sep="")

  }

  if (!is.null(expect.cells)){
   params<-paste(params," --expect-cell=",expect.cells, sep="")
  }

  if (!is.null(force.cells)){
   params<-paste(params," --force-cells=",force.cells, sep="")
  }

  if (nosecondary){
   params<-paste(params," --nosecondary", sep="")
  }

  if (!is.null(chemistry)){
   params<-paste(params," --chemistry=",chemistry, sep="")
  }

  if (!is.null(r1.length)){
   params<-paste(params," --r1-length=",r1.length, sep="")
  }

  if (!is.null(r2.length)){
   params<-paste(params," --r2-length=",r2.length, sep="")
  }

  if (!is.null(lanes)){
   params<-paste(params," --lanes=",lanes, sep="")
  }

  if (!is.null(localcores)){
   params<-paste(params," --localcores=",localcores, sep="")
  }

  if (!is.null(localmem)){
   params<-paste(params," --localmem=", localmem, sep="")
  }

  params.split <- unlist(strsplit(params, dockerImage))
  params0 <- paste(params.split[1], " ", dockerImage, " /bin/sh /data/script.sh", sep="")
  cat(params0,"\n")
  params1 <- NULL
  params1[1] <- "cd /data"
  params1[2] <- params.split[2]
  params1[3] <- paste("chmod 777 -R /data/", id, sep="")
  if(version=="2"){
    params1[4] <- paste("/bin/cellranger mat2csv /data/", id,"/outs/filtered_gene_bc_matrices ",id,".csv", sep="")
  }else if(version=="3"){
    params1[4] <- paste("/bin/cellranger mat2csv /data/", id,"/outs/filtered_feature_bc_matrix ",id,".csv", sep="")
  }
  



  fileConn<-file(paste(scrat_tmp.folder,"/script.sh", sep=""), "w")
  writeLines(params1, fileConn)
  close(fileConn)
  system(paste("chmod +x ", scrat_tmp.folder,"/script.sh", sep=""))


  #Run docker
  resultRun <- runDocker(group=group, params=params0)
  #waiting for the end of the container work
  if(resultRun==0){
    system(paste("cp -R ", scrat_tmp.folder, "/", id, " ", fastq.folder, sep=""))
    system(paste("cp ", scrat_tmp.folder, "/results_cellranger.csv ", fastq.folder, sep=""))
    system(paste("sed \'s|,|\t|g\' ",fastq.folder,"/",id,".csv > ", fastq.folder,"/",id,".txt", sep=""))
    cat("\nCellranger analysis is finished\n")
  }

 #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(fastq.folder)
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
  container.id <- readLines(paste(fastq.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", substr(container.id,1,12), " &> ",fastq.folder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  #removing temporary folder
#  cat("\n\nRemoving the temporary file ....\n")
#  system(paste("rm -fR ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",fastq.folder, sep=""))
  setwd(home)

}

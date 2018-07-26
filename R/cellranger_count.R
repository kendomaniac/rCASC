#' @title Cellranger count
#' @description This function takes FASTQ files from cellranger mkfastq and performs alignment, filtering, barcode counting, and UMI counting. It uses the Chromium cellular barcodes to generate gene-barcode matrices, determine clusters, and perform gene expression analysis. The count pipeline can take input from multiple sequencing runs on the same library.

#' --id	     A unique run ID string: e.g. sample345
#' --fastq	 Either: Path of the fastq_path folder generated by cellranger mkfastq
#'           e.g. /home/jdoe/runs/HAWT7ADXX/outs/fastq_path. This contains a directory hierarchy that cellranger count will              #'           automatically traverse.
#'           - OR -
#'           Any folder containing fastq files, for example if the fastq files were generated by a service provider and delivered   #'           outside the context of the mkfastq output directory structure.
#'           Can take multiple comma-separated paths, which is helpful if the same library was sequenced on multiple flowcells.
#'           Doing this will treat all reads from the library, across flowcells, as one sample.
#'           If you have multiple libraries for the sample, you will need to run cellranger count on them individually, and then #'           combine them with cellranger aggr.
#'--sample	 Sample name as specified in the sample sheet supplied to cellranger mkfastq.
#'           Can take multiple comma-separated values, which is helpful if the same library was sequenced on multiple flowcells #'           and the sample name used (and therefore fastq file prefix) is not identical between them.
#'           Doing this will treat all reads from the library, across flowcells, as one sample.
#'           If you have multiple libraries for the sample, you will need to run cellranger count on them individually, and then #'           combine them with cellranger aggr.
#'           Allowable characters in sample names are letters, numbers, hyphens, and underscores.
#'--transcriptome       	Path to the Cell Ranger compatible transcriptome reference e.g.
#'                          For a human-only sample, use /opt/refdata-cellranger-GRCh38-1.2.0
#'                          For a human and mouse mixture sample, use /opt/refdata-cellranger-hg19-and-mm10-1.2.0
#'--expect-cells    (optional) Expected number of recovered cells. Default: 3,000 cells.
#'--force-cells	    (optional) Force pipeline to use this number of cells, bypassing the cell detection algorithm.
#'                  Use this if the number of cells estimated by Cell Ranger is not consistent with the barcode rank plot.
#'--nosecondary	    (optional) Add this flag to skip secondary analysis of the gene-barcode matrix (dimensionality reduction,
#'                  clustering and visualization). Set this if you plan to use cellranger reanalyze or your own custom analysis.
#'--chemistry	    (optional) Assay configuration. One of:
#'                  auto for autodetection (default), threeprime for Single Cell 3′, fiveprime for Single Cell 5′, SC3Pv1 for
#'                  Single Cell 3′ v1, SC3Pv2 for Single Cell 3′ v2, SC5P-PE for Single Cell 5′ paired-end (both R1 and R2 are
#'                  used for alignment),SC5P-R2 for Single Cell 5′ R2-only (where only R2 is used for alignment).
#'--r1-length	    (optional) Hard-trim the input R1 sequence to this length. Note that the length includes the Barcode and UMI
#'                  sequences so do not set this below 26 for Single Cell 3′ v2 or Single Cell 5′. This and --r2-length are
#'                  useful for determining the optimal read length for sequencing.
#'--r2-length	    (optional) Hard-trim the input R2 sequence to this length.
#'--lanes	        (optional) Lanes associated with this sample
#'--localcores	    Restricts cellranger to use specified number of cores to execute pipeline stages. By default, cellranger
#'                  will use all of the cores available on your system.
#'--localmem	    Restricts cellranger to use specified amount of memory (in GB) to execute pipeline stages. By default,
#'                  cellranger will use 90% of the memory available on your system. Please note that cellranger requires at least
#'                  16 GB of memory to run all pipeline stages.
#'--indices	        (Deprecated. Optional. Only used for output from cellranger demux) Sample indices associated with this
#'                  sample. Comma-separated list of:
#'                  index set plate well: SI-3A-A1
#'                  index sequences: TCGCCATA,GTATACAC



#' @author ????
#'
#' @examples
#' \dontrun{
#' library(casc)
#' downloadContainers(group="docker","?????")
#' system("wget ??????")
#' AGGIUNGERE ESEMPIO CON POSSIBILITA' DI SCARICARE DATO
#'
#'
#'
#'
#' }
#'
#'
#' @export

cellranger_count <- function(group=c("sudo","docker"),  id, transcriptome.folder,  fastq.folder,  sample, expect.cells=NULL, force.cells=NULL, nosecondary=NULL, chemistry=NULL, r1.length=NULL,  r2.length=NULL, lanes=NULL, localcores=NULL, localmem=NULL,  indices=NULL, scratch.folder){

  #docker image
  dockerImage="docker.io/grromano/cellranger"

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







  #executing the docker job
  #Le directory vanno montate tutte con il -v  user:doker
  #modifica qui /bin/checkscript.sh
  params <- paste("--cidfile ",fastq.folder,"/dockerID -v ",transcriptome.folder,":/transcr -v ", fastq.folder, ":/data -d ",dockerImage, " /bin/cellranger count  --id=",id," --transcriptome=/transcr --fastqs=/data", sep="")

  cat(params,"\n")

  if(!is.null(sample)){

    params<-paste(params,"--sample",sample)

  }

  if (!is.null(expect.cells)){
   params<-paste(params," --expect-cell=",expect.cells)
  }

  if (!is.null(force.cells)){
   params<-paste(params," --force-cells=",force.cells)
  }

  if (!is.null(nosecondary)){
   params<-paste(params," --nosecondary=",nosecondary)
  }

  if (!is.null(chemistry)){
   params<-paste(params," --chemistry=",chemistry)
  }

  if (!is.null(r1.length)){
   params<-paste(params," --r1-length=",r1.length)
  }

  if (!is.null(r2.length)){
   params<-paste(params," --r2-length=",r2.length)
  }

  if (!is.null(lanes)){
   params<-paste(params," --lanes=",lanes)
  }

  if (!is.null(localcores)){
   params<-paste(params," --localcores=",localcores)
  }

  if (!is.null(localmem)){
   params<-paste(params," --localmem=", localmem)
  }

  if (!is.null(indices)){
   params<-paste(params," --indices=",indices)
  }


  #Run docker
  resultRun <- runDocker(group=group, params=params)
  #waiting for the end of the container work
  if(resultRun==0){
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
  container.id <- readLines(paste(fastq.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", substr(container.id,1,12), " &> ",fastq.folder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -fR ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="casc"),"containers/containers.txt",sep="/")," ",fastq.folder, sep=""))
  setwd(home)

}
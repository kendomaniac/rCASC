#' @title A function to create a genome index for indrop V2 single cell data
#' @description This function executes a docker that produces as output the genome index index for bowtie
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param index.folder, a character string indicating the folder where the index will be created. The index will have the prefix genome.
#' @param ensembl.urlgenome, a character string indicating the URL from ENSEMBL ftp for the unmasked genome sequence of interest
#' @param ensembl.urlgtf, a character string indicating the URL from ENSEMBL ftp for the GTF for genome of interest
#' @author Raffaele Calogero and Riccardo Panero, raffaele.calogero [at] unito [dot] it, Bioinformatics and Genomics unit, University of Torino Italy
#'
#' @examples
#' \dontrun{

#' library(casc)
#' #running indropCounts index build
#' indropIndex(group="docker", index.folder=getwd(),
#'     ensembl.urlgenome="ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz",
#'     ensembl.urlgtf="ftp://ftp.ensembl.org/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf.gz")
#' }
#'
#' @export
indropIndex <- function(group=c("sudo","docker"), index.folder, ensembl.urlgenome, ensembl.urlgtf){

  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    return()
  }
  #storing the position of the home folder
  home <- getwd()
  #running time 1
  ptm <- proc.time()
  #setting the data.folder as working folder
  if (!file.exists(index.folder)){
    cat(paste("\nIt seems that the ",index.folder, " folder does not exist\n"))
    return(2)
  }
  setwd(index.folder)



  yaml.file=paste(path.package(package="casc"),"extras/indrop.yaml",sep="/")
  system(paste("cp ",yaml.file," ", file.path(index.folder),sep=""))
  system(paste("chmod 777 -R", file.path(index.folder)))

  #edit yaml
  yaml <- readLines("indrop.yaml")
  project_name <- yaml[grep("project_name", yaml)]
  project_name <- sub("CRISPR", "INDEX", project_name)
  yaml[grep("project_name", yaml)] <- project_name

#  project_dir <- yaml[grep("project_dir", yaml)]
#  project_dir <- sub("/sto2/labcamargo/Documents/single_cell/CRISPR_single_cell_9Nov17/inDrops/", "/data/scratch", project_dir)
#  yaml[grep("project_dir", yaml)] <- project_dir

#  sample_name <- yaml[grep("  - name :", yaml)]
#  sample_name <- sub("CRISPR", sample.name, sample_name)
#  yaml[grep("  - name :", yaml)] <- sample_name

#  input_dir <- yaml[grep("    dir :", yaml)]
#  input_dir <- sub("/sto2/labcamargo/Documents/single_cell/CRISPR_single_cell_9Nov17/basespace/171004_M00620_0217_000000000-BFWPC_FASTQ", "/data/scratch/input", input_dir)
#  yaml[grep("    dir :", yaml)] <- input_dir

#  split_affixes <- yaml[grep("    split_affixes :", yaml)]
#  split_affixes <- sub("S1_L001", split.affixes, split_affixes)
#  yaml[grep("    split_affixes :", yaml)] <- split_affixes

#  library_name <- yaml[grep("library_name:", yaml)]
#  library_name <- gsub("Sample1", sample.name, library_name)
#  yaml[grep("library_name:", yaml)] <- library_name

  bowtie_index <- yaml[grep("bowtie_index :", yaml)]
  bowtie_index <- gsub("/sto2/labcamargo/Documents/bowtie_index/mm10/Mus_musculus.GRCm38.85.index", "/index/genome", bowtie_index)
  yaml[grep("bowtie_index :", yaml)] <- bowtie_index

  #UMI parameters
#  m <- yaml[grep("    m : 10 #Ignore reads with more than M alignments, after filtering on distance from transcript end.", yaml)]
#  m <- sub("10", M, m)
#  yaml[grep("    m : 10 #Ignore reads with more than M alignments, after filtering on distance from transcript end.", yaml)] <- m

#  u <- yaml[grep("    u : 2 #Ignore counts from UMI that should be split among more than U genes.", yaml)]
#  u <- sub("2", U, u)
#  yaml[grep("    u : 2 #Ignore counts from UMI that should be split among more than U genes.", yaml)] <- u

#  d <- yaml[grep("    d : 400 #Maximal distance from transcript end, NOT INCLUDING THE POLYA TAIL", yaml)]
#  d <- sub("400", D, d)
#  yaml[grep("    d : 400 #Maximal distance from transcript end, NOT INCLUDING THE POLYA TAIL", yaml)] <- d

  #outout params
#  low_complexity_mask <- yaml[grep("    low_complexity_mask: False", yaml)]
#  low_complexity_mask <- sub("False", low.complexity.mask, low_complexity_mask)
#  yaml[grep("    low_complexity_mask: False", yaml)] <- low_complexity_mask


  zz <- file("indrop.yaml", "w")
  writeLines(yaml, zz)
  close(zz)
  #

  system(paste("chmod 777 -R", file.path(index.folder)))

  params <- paste("--cidfile ",index.folder,"/dockerID -v ", index.folder,":/index -d docker.io/repbioinfo/indrop.2017.01 sh /bin/indropIndex.sh ", ensembl.urlgenome, " ", ensembl.urlgtf, sep="")
  resultRun <- runDocker(group=group, params=params)

#  if(group=="sudo"){
#    params <- paste("--cidfile ",index.folder,"/dockerID -v ", index.folder,":/index -d docker.io/repbioinfo/indrop.2017.01 sh /bin/indropIndex.sh ", ensembl.urlgenome, " ", ensembl.urlgtf, sep="")
#    resultRun <- runDocker(group="sudo",container="docker.io/repbioinfo/indrop.2017.01", params=params)
#  }else{
#    params <- paste("--cidfile ",index.folder,"/dockerID -v ", index.folder,":/index -d docker.io/repbioinfo/indrop.2017.01 sh /bin/indropIndex.sh ", ensembl.urlgenome, " ", ensembl.urlgtf, sep="")
#    resultRun <- runDocker(group="docker",container="docker.io/repbioinfo/indrop.2017.01", params=params)
#  }


  if(resultRun==0){
    cat("\n inDrop genome index generation is finished\n")
  }

  #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(index.folder)
  dir <- dir[grep("run.info",dir)]
  if(length(dir)>0){
    con <- file("run.info", "r")
    tmp.run <- readLines(con)
    close(con)
    tmp.run[length(tmp.run)+1] <- paste("inDrop index user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("inDrop index system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("inDrop index elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("inDrop index user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("inDrop index system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("inDrop index elapsed run time mins ",ptm[3]/60, sep="")

    writeLines(tmp.run,"run.info")
  }

  #saving log and removing docker container
  container.id <- readLines(paste(index.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", substr(container.id,1,12), " &> ",index.folder,"/indropIndex_", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))

  cat("\n\nRemoving the temporary file ....\n")
  system("rm -fR dockerID")
  system(paste("cp ",paste(path.package(package="casc"),"containers/containers.txt",sep="/")," ",index.folder, sep=""))
  setwd(home)

}

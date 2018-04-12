#' @title A function to handle sigle cell Lorenz Quality filter for Single-cells
#' @description This function executes a docker that produces ....
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param data.folder, a character string indicating the folder where input data are located and where output will be written
#' @param matrixName, counts table name. Matrix data file must be in data.folder. The file MUST contain RAW counts, without any modification, such as log transformation, normalizatio etc.
#' @param umiXgene, a integer defining how many UMI are required to call a gene as present. default: 3
#' @param p_value, threshold to be used for the filtering
#' @param format, counts file extension, "txt", "csv"
#' @param separator, separator used in count file, e.g. '\\t', ','
#'
#' @author Name Family name, myemail [at] somewhere [dot] org, Affiliation
#'
#' @return output will be in the same format and with the same separator of input. It also returns a PDF with the genes vs UMI plot, discarded_cells.pdf. In blue discarded cells are shown.
#'
#' @examples
#' \dontrun{
#'  system("wget http://130.192.119.59/public/lorenz.zip")
#'  unzip("lorenz.zip")
#'  setwd("./lorenz")
#'  library("CASC")
#'  lorenzFilter(group="docker", scratch.folder="/data/scratch",
#'           data.folder=getwd(), matrixName="Buettner", p_value=0.05, umiXgene=3, format="csv", separator=',')
#' }
#'
#' @export
lorenzFilter <- function(group=c("sudo","docker"), scratch.folder, data.folder, matrixName, p_value, umiXgene=3, format, separator){
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
  if (!file.exists(data.folder)){
    cat(paste("\nIt seems that the ",data.folder, " folder does not exist\n"))
    return(2)
  }
  setwd(data.folder)
  #check  if scratch folder exist
  if (!file.exists(scratch.folder)){
    cat(paste("\nIt seems that the ",scratch.folder, " folder does not exist\n"))
    return(3)
  }
  tmp.folder <- gsub(":","-",gsub(" ","-",date()))
  scrat_tmp.folder=file.path(scratch.folder, tmp.folder)
  writeLines(scrat_tmp.folder,paste(data.folder,"/tempFolderID", sep=""))
  cat("\ncreating a folder in scratch folder\n")
  dir.create(file.path(scrat_tmp.folder))

  write.csv(read.table(paste(data.folder,"/",matrixName,".",format,sep=""),sep=separator,header=TRUE,row.names=1),paste(scrat_tmp.folder,"/set1.csv",sep=""))
if(separator=="\t"){
separator="tab"
}

  #executing the docker job
  if(group=="sudo"){
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ", data.folder, ":/data -d docker.io/rcaloger/lorenz Rscript /home/main.R ",matrixName," ",p_value," ",format," ",separator, sep="")
    resultRun <- runDocker(group="sudo",container="docker.io/rcaloger/lorenz", params=params)
  }else{
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ", data.folder, ":/data -d docker.io/rcaloger/lorenz Rscript /home/main.R ",matrixName," ",p_value," ",format," ",separator, sep="")
    resultRun <- runDocker(group="docker",container="docker.io/rcaloger/lorenz", params=params)
  }
  #waiting for the end of the container work
  if(resultRun==0){
    system(paste("cp ", scrat_tmp.folder, "/* ", data.folder, sep=""))
  }
  if(separator=="tab"){
    separator="\t"
  }
  dir <- dir(data.folder)
  files <- dir[grep(matrixName, dir)]
  if(length(files) == 2){
      files.lorenz <- dir[grep("^lorenz", dir)]
      output <- intersect(files, files.lorenz)
      input <- setdiff(files, files.lorenz)
      #plotting the genes vs umi all cells
      tmp0 <- read.table(input, sep=separator, header=T, row.names=1)
      genes <- list()
      for(i in 1:dim(tmp0)[2]){
        x = rep(0, dim(tmp0)[1])
        x[which(tmp0[,i] >=  umiXgene)] <- 1
        genes[[i]] <- x
      }
      genes <- as.data.frame(genes)
      genes.sum <-  apply(genes,2, sum)
      umi.sum <- apply(tmp0,2, sum)
      pdf("discarded-cells.pdf")
      plot(log10(umi.sum), genes.sum, xlab="log10 UMI", ylab="# of genes")
      points(log10(umi.sum), genes.sum, pch=19, cex=0.5, col="blue")

      tmp <- read.table(output, sep=separator, header=T, row.names=1)
      genes <- list()
      for(i in 1:dim(tmp)[2]){
        x = rep(0, dim(tmp)[1])
        x[which(tmp[,i] >=  umiXgene)] <- 1
        genes[[i]] <- x
      }
      genes <- as.data.frame(genes)
      genes.sum <-  apply(genes,2, sum)
      umi.sum <- apply(tmp,2, sum)
      points(log10(umi.sum), genes.sum, pch=19, cex=0.5, col="red")
      dev.off()

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
  system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/lorenzFilter_", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="casc"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}

#' @title A function to convert raw count in log counts
#' @description This function convert raw counts in log format
#' @param file, a character string indicating the path of the file. IMPORTANT: full path to the file MUST be included.
#' @param log.base, the base of the log to be used for the transformation
#' @param type, the type of input file, txt, tab delimited. csv, comma separated
#' @author Raffaele Calogero, raffaele.calogero [at] unito [dot] it, University of Torino
#' @return log transformed table
#' @examples
#' \dontrun{
#'     system("wget http://130.192.119.59/public/TO BE INSERTED")
#'     #running skeleton
#'     counts2log(group="docker", counts.matrix, data.folder=getwd(), log.base=10)
#' }
#'
#' @export

counts2log <- function(file, log.base=c(2,10)){

  home <- getwd()
  data.folder=dirname(file)
  setwd(data.folder)
  positions=length(strsplit(basename(file),"\\.")[[1]])
  matrixNameC=strsplit(basename(file),"\\.")[[1]]
  matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
  format=strsplit(basename(basename(file)),"\\.")[[1]][positions]
  counts.matrix <- paste(matrixName, format, sep=".")

  if(format=="txt"){
    tmp <- read.table(counts.matrix, sep="\t", stringsAsFactors = FALSE, header=T, check.names = FALSE, row.names=1)
    if(log.base==2){
      tmpl <- log2(tmp +1)
      write.table(tmpl, paste("log2_",sub(".txt","", counts.matrix), ".csv", sep=""), sep=",", col.names=NA)
    }else if(log.base==10){
      tmpl <- log10(tmp +1)
      write.table(tmpl, paste("log10_",sub(".txt","", counts.matrix), ".csv", sep=""), sep=",", col.names=NA)
    }
  }else if(format=="csv"){
    tmp <- read.table(counts.matrix, sep=",", stringsAsFactors = TRUE, header=T, check.names = FALSE, row.names=1)
    if(log.base==2){
      tmpl <- log2(tmp +1)
      write.table(tmpl, paste("log2_",sub(".csv","", counts.matrix), ".csv", sep=""), sep=",", col.names=NA)
    }else if(log.base==10){
      tmpl <- log10(tmp +1)
      write.table(tmpl, paste("log10_",sub(".csv","", counts.matrix), ".csv", sep=""), sep=",", col.names=NA)
    }
  }
  setwd(home)
}

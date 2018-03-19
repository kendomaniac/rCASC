#' @title A function to convert raw count in log counts
#' @description This function convert raw counts in log format
#' @param counts.matrix, a character string indicating the name of the input raw counts file.
#' @param data.folder, a character string indicating the folder where input data are located and where output will be written
#' @param log.base, the base of the log to be used for the transformation
#' @param type, the type of input file, txt, tab delimited. csv, comma separated
#' @author Raffaele Calogero, raffaele.calogero [at] unito [dot] it, University of Torino
#' @return log transformed table
#' @examples
#' \dontrun{
#'     system("wget http://130.192.119.59/public/test_R1.fastq.gz")
#'     #running skeleton
#'     fastqc(group="docker", data.folder=getwd())
#' }
#'
#' @export

counts2log <- function(counts.matrix, data.folder=getwd(), log.base=c(2,10), type=c("txt","csv")){
  if(type=="txt"){
    tmp <- read.table(counts.matrix, sep="\t", stringsAsFactors = FALSE, header=T, check.names = FALSE, row.names=1)
    if(log.base==2){
      tmpl <- log2(tmp +1)
      write.table(tmpl, paste("log2_",sub(".txt","", counts.matrix), ".csv", sep=""), sep=",", col.names=NA)
    }else if(log.base==10){
      tmpl <- log10(tmp +1)
      write.table(tmpl, paste("log10_",sub(".txt","", counts.matrix), ".csv", sep=""), sep=",", col.names=NA)
    }
  }else if(type=="csv"){
    tmp <- read.table(counts.matrix, sep=",", stringsAsFactors = TRUE, header=T, check.names = FALSE, row.names=1)
    if(log.base==2){
      tmpl <- log2(tmp +1)
      write.table(tmpl, paste("log2_",sub(".csv","", counts.matrix), ".csv", sep=""), sep=",", col.names=NA)
    }else if(log.base==10){
      tmpl <- log10(tmp +1)
      write.table(tmpl, paste("log10_",sub(".csv","", counts.matrix), ".csv", sep=""), sep=",", col.names=NA)
    }
  }
}

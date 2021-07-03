% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixcr.R
\name{mixcr}
\alias{mixcr}
\title{mixcr}
\usage{
mixcr(
  group = c("sudo", "docker"),
  scratch.folder,
  fastqPath,
  resFolderCustom = "NULL"
)
}
\arguments{
\item{group, }{a character string. Two options: sudo or docker, depending to which group the user belongs}

\item{scratch.folder, }{a character string indicating the path of the scratch folder}

\item{fastqPath, }{Path of fastq folder}

\item{resFolderCustom, }{optional parameter. Default will store the results in fastqPath otherwise will store the results in resFolderCustom path.}
}
\value{
an indexed genome compliant with 10XGenomics cellranger
}
\description{
This function creates mixcrFiles
}
\examples{
\dontrun{
library(rCASC)
dir.create("scratch")
scratch.folder=paste(getwd(),"scratch",sep="/")
fastqPath=paste(getwd(),"fastq",sep="/")
resFolder=paste(getwd(),"resFolder",sep="/")
dir.create(resFolder)
mixcr(group=c("sudo"),scratch.folder=scratch.folder,fastqPath=fastqPath,resFolderCustom=resFolder)
}


}
\author{
Luca Alessandrì
}
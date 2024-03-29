% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autoanalysis4Clustering.R
\name{autoencoderAnalysis4Clustering}
\alias{autoencoderAnalysis4Clustering}
\title{Autoencoder Analysis}
\usage{
autoencoderAnalysis4Clustering(
  group = c("sudo", "docker"),
  scratch.folder,
  file,
  separator,
  projectName,
  seed = 1111,
  Sp
)
}
\arguments{
\item{group, }{a character string. Two options: sudo or docker, depending to which group the user belongs}

\item{scratch.folder, }{a character string indicating the path of the scratch folder}

\item{file, }{a character string indicating the path of the file, with file name and extension included.Has to be the one in the projectName folder.Different so from the previous one.}

\item{separator, }{separator used in count file, e.g. '\\t', ','}

\item{projectName, }{might be different from the matrixname in order to perform different analysis on the same dataset}

\item{seed, }{important value to reproduce the same results with same input}

\item{Sp, }{minimun number of percentage of cells that has to be in common between two permutation to be the same cluster.}
}
\description{
This function Compress data using autoencoder partially connected
}
\examples{
\dontrun{
 autoencoderAnalysis4Clustering(group=group, scratch.folder=scratch.folder, file=file2,separator=separator,projectName=projectName2,seed=seed,Sp=Sp)
}
}
\author{
Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
}

#' @title Download for the first time all containers embedded in the workflows
#' @description This is a functin that preapre the docker environment to be used for the first time the docker4seq is installed.
#' @param group, a character string. Two options: \code{"sudo"} or \code{"docker"}, depending to which group the user belongs
#' @param containers.file, a character string with the name of the file which indicate which are the initial set of containers to be downloaded. options: full, mini, nano, sca The set is given by a file located in the folder containers of docker4seq package, full indicates a complete installation, mini refers to an installation including few preprocessing and all clustering tools, nano even a smaller implementation, sca the minimal subset of dockers to run an analysis with Partialy Connected Autoencoders
#' @author Raffaele Calogero
#'
#' @examples
#'\dontrun{
##'     #running runDocker
#'      downloadContainers(group="docker", containers.file="full")
#'
#' }
#' @export
downloadContainers <- function(group="docker", containers.file=c("full", "mini", "nano", "sca"){
   if(containers.file=="full"){
     containers.file=paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")
   }else if(containers.file=="mini"){
     containers.file=paste(path.package(package="rCASC"),"containers/mini.txt",sep="/")
   }
   else if(containers.file=="nano"){
      containers.file=paste(path.package(package="rCASC"),"containers/nano.txt",sep="/")
   }
   else if(containers.file=="sca"){
      containers.file=paste(path.package(package="rCASC"),"containers/sca_only.txt",sep="/")
   }
   containers <- readLines(containers.file)
   for(i in containers){
     if(group=="sudo"){
       system(paste("sudo docker pull ",i, sep=""))
     }else if(group=="docker"){
       system(paste("docker pull ",i, sep=""))
     }else{
       cat("\nThe group provides is neither sudo or docker\n")
       return(2)
     }
   }
   writeLines(containers, paste(path.package(package="rCASC"),"containers/containers.txt",sep="/"))
}

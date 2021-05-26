#' @title gseaXLmHG
#' @description This function execute XLmHG estimation of the presence of enriched GSEA class in cluster-specific genes detected by COMET for the clusters of two independent dataset and uses this information to search for correspondence between clusters in the two independent datasets. This function requires that the two datasets have been clustered in rCASC and COME analysos was run.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param xCometFolder, path of Comet results from X experiment
#' @param yCometFolder, path of Comet results from Y experiment
#' @param gsea, list of the available GSEA classes: c1.all, c2.cgp, c2.cp.biocarta, c2.cp.kegg, c2.cp.pid, c2.cp.reactome, c2.cp.wikipathways, c3.all, c3.mir, c3.tft.gtrd, c3.tft, c4.cgn, c4.cm, c5.go.bp, c5.go.cc, c5.go.mf, c5.hpo, c6.all, c7.all, c8.all, h.all, msigdb.all. Please note that msigdb.all includes all gsea classes.
#' @param X, X parameter for the XLmHG, default 5, for more info please see XLmHG help: https://xl-mhg.readthedocs.io/en/latest/.
#' @param L, L parameter for the XLmHG, default 0.15, for more info please see XLmHG help: https://xl-mhg.readthedocs.io/en/latest/.
#' @param pvalue, XLmHG pvalue threshold
#' @param separatorX, separator used in count file, e.g. '\\t', ','
#' @param separatorY, separator used in count file, e.g. '\\t', ','
#' @param outputFolder, the folder where GSEA folder will be created
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return A folder called GSEA with all the results generated. The function produces an integrated output, final_score.csv, of the results obtained for progressively increasing length of the cluster-specific-genes detected by COMET. 
#' @examples
#' \dontrun{
#'  library(rCASC)
#'  gseaXLmHG(group="docker",
#'         scratch.folder="/scratch", 
#'         xCometFolder="/data/reanalysis_on_AIsc/comparing_CRC0327/NT_CTX/CRC0327_NT_2_clx/VandE/Results/VandE/8/outputdata",
#'         yCometFolder="/data/reanalysis_on_AIsc/comparing_CRC0327/NT_CTX/CRC0327_cetux_2_clx/VandE/Results/VandE/8/outputdata",
#'         gsea="c2.cp.kegg",
#'         X=5,
#'         L=0.15,
#'         pvalue=0.05,
#'         separatorX=",",
#'         separatorY=",",
#'         outputFolder="/data/reanalysis_on_AIsc/comparing_CRC0327/NT_CTX"
#'  )
#'}
#' @export
gseaXLmHG <- function(group=c("sudo","docker"), scratch.folder, xCometFolder, yCometFolder, 
                   gsea=c("c1.all","c2.cgp","c2.cp.biocarta","c2.cp.kegg","c2.cp.pid",
                          "c2.cp.reactome","c2.cp.wikipathways", "c3.all","c3.mir","c3.tft.gtrd",
                          "c3.tft","c4.cgn","c4.cm","c5.go.bp","c5.go.cc","c5.go.mf","c5.hpo",
                          "c6.all","c7.all","c8.all","h.all","msigdb.all"), X, L, pvalue, separatorX, separatorY, outputFolder){
 
  data.folder=outputFolder
  #running time 1
  ptm <- proc.time()
  #setting the data.folder as working folder
  if (!file.exists(data.folder)){
    cat(paste("\nIt seems that the ",data.folder, " folder does not exist\n"))
    return(2)
  }

  #storing the position of the home folder
  home <- getwd()
  setwd(data.folder)
  #initialize status
  system("echo 0 > ExitStatusFile 2>&1")

  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    system("echo 10 > ExitStatusFile 2>&1")
    setwd(home)
    return(10)
  }



  #check  if scratch folder exist
  if (!file.exists(scratch.folder)){
    cat(paste("\nIt seems that the ",scratch.folder, " folder does not exist\n"))
    system("echo 3 > ExitStatusFile 2>&1")
    setwd(data.folder)
    return(3)
  }
  tmp.folder <- gsub(":","-",gsub(" ","-",date()))
  scrat_tmp.folder=file.path(scratch.folder, tmp.folder)
  writeLines(scrat_tmp.folder,paste(data.folder,"/tempFolderID", sep=""))
  cat("\ncreating a folder in scratch folder\n")
  dir.create(file.path(scrat_tmp.folder))
  #preprocess matrix and copying files

if(separatorX=="\t"){
separatorX="tab"
}

if(separatorY=="\t"){
separatorY="tab"
}
  system(paste("cp -r ",xCometFolder," ",scrat_tmp.folder,"/",sep=""))
  system(paste("mv ", scrat_tmp.folder,"/outputdata ", scrat_tmp.folder,"/Xoutputdata", sep=""))

  system(paste("cp -r ",yCometFolder," ",scrat_tmp.folder,"/",sep=""))
  system(paste("mv ", scrat_tmp.folder,"/outputdata ", scrat_tmp.folder,"/Youtputdata", sep=""))

  xCometFolder="/data/Xoutputdata"
  yCometFolder="/data/Youtputdata"

  #executing the docker job

  params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/data -d docker.io/repbioinfo/xlmhg.2021.01 Rscript /home/loop_toprnk.R ",gsea," ", X, " ", L, " ", pvalue, " ",xCometFolder," ",yCometFolder," ",separatorX," ",separatorY, sep="")

resultRun <- runDocker(group=group, params=params)

  #waiting for the end of the container work
  if(resultRun==0){
    #system(paste("cp ", scrat_tmp.folder, "/*top* ", data.folder, sep=""))
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
  system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))


  #Copy result folder
  cat("Copying Result Folder")
    system(paste("cp -r ",scrat_tmp.folder,"/GSEA ",data.folder,"/",sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}
 

 
#' @title BCscTool_wrapper															
#' @description The present function associates the clusters of two different and indipendent experiments using the Bray Curtis dissimilarity. Starting from the two original sets it takes only the genes common to both and calculates the Bray Curtis Dissimilarity between the clusters of one set and those of the other. all the results will be saved in the folder called "BCsc" inside the folder where the files of the first set are located 	
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs					
#' @param setX, path to the first dataset
#' @param nclustX, number of cluster of the first dataset.
#' @param markerX, path to the FOLDER containing cluster-specific genes for the FIRST dataset. (after using Comet directory '.../nclustX/outputdata') 
#' @param setY, path to the second dataset
#' @param nclustY, number of cluster of the second dataset
#' @param markerY, path to the FOLDER containing cluster-specific genes for the SECOND dataset. (after using Comet directory '.../nclustY/outputdata')
#' @param separator, separator used in both count file, e.g. '\\t', ','
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param output.folder, a character string indicating the path of the output folder to save analysis
#' @param pvalueComet, parameter to select only marker genes  with a pvalue under a threshold value. DEFAULT = 1
#' @param starting_genes, number of starting genes for each cluster selected from the main marker gene lists. DEFAULT = 500
#' @param subdivision, number of progressive reductions of the starting genes list. DEFAULT = 10 
#' @param contamination, percentage of genes to be contaminated at each permutation during the Bray-Curtis calculation. DEFAULT = 5
#' @param permutation, number of Bray-Curtis permutations. DEFAULT = 50
#' @param threshold, minimum value to consider a gene expressed in both count sets. DEFAULT VALUE: 3
#' @author Gabriele Piacenti, g [dot] pia91 [at] gmail [dot] com, University of Torino
#'
#' @return 
#' @examples			
#' \dontrun{
#' BCscWrapper(group="docker",
#'       setX = '/2tb/torino/piacenti/NT1_NT2/CRC0327_NT_1_clx/VandE/VandE.csv',
#'       nclustX = 7,
#'       markerX = '/2tb/torino/piacenti/NT1_NT2/CRC0327_NT_1_clx/VandE/Results/VandE/7/outputdata',
#'       setY = '/2tb/torino/piacenti/NT1_NT2/CRC0327_NT_2_clx/VandE/VandE.csv',
#'       nclustY = 8,
#'      markerY = '/2tb/torino/piacenti/NT1_NT2/CRC0327_NT_2_clx/VandE/Results/VandE/8/outputdata',
#'      scratch.folder = '/2tb/torino/piacenti/user09/piacenti/Prova_Prove/Script/scratch',
#'      output.folder = '/2tb/torino/piacenti/user09/piacenti/Prova_Prove/Script',
#'      pvalueComet = 1,
#'      subdivision = 10,
#'      starting_genes = 500,
#'      contamination = 5,
#'      permutation = 50,
#'      threshold = 3)
#' }
#' @export
BCscWrapper <- function(group=c("sudo","docker"),
                        setX,
                        nclustX,
                        markerX,
                        setY,
                        nclustY,
                        markerY,
                        separator = ",",
                        scratch.folder,
                        output.folder, 
                        pvalueComet = 1,
                        subdivision = 10,
                        starting_genes = 500,
                        contamination = 5,
                        permutation = 50,
                        threshold = 3) {   
                        

data.folder=output.folder
dir.create(output.folder,showWarnings = F)

positions=length(strsplit(basename(setX),"\\.")[[1]])
matrixNameC=strsplit(basename(setX),"\\.")[[1]]
matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
format=strsplit(basename(basename(setX)),"\\.")[[1]][positions]


## da qui fino a segno *** NON TOCCARE

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

if(separator=="\t"){
separator="tab"
}

#### ***


Marker <- function(clustering, path, pvalue, scratch, separator){
  
  dir <- scratch
  name <- strsplit(basename(clustering),".csv")[[1]][1]
  
  cluster_SETA <- read.csv(file = paste(scratch,clustering,sep = "/"),
                           header = TRUE,
                           sep = separator,row.names = 1)
  
  cluster_SETA <- sort(unique(cluster_SETA$Belonging_Cluster))
  
  gene <- c()
  cluster <-c()
  pval <-c()
  path <- path
  
  for(n in 1:length(cluster_SETA)){
    
    table <- read.csv(file = paste(path,
                                   "/cluster_",
                                   cluster_SETA[n],
                                   "_singleton_positive_markers_ranked.csv",
                                   sep="" ),
                      header = TRUE,
                      sep = separator,
                      row.names = 1)
    
    table <- table[table$q_value <= pvalue,]
    
    gene <- c(gene, rownames(table))
    cluster <- c(cluster, rep(cluster_SETA[n], nrow(table)))
    pval <- c(pval, table$q_value)
    
    
  }
  
  markergenes <- data.frame(gene,cluster,pval)
  
  write.table(markergenes,
              file = paste(dir,"/",name,"_MINOR_",pvalue,".marker.csv",sep =""),
              sep = separator,
              col.names = NA)
  
}

set_x <- read.table(file =setX, header = TRUE,sep = separator,row.names = 1)
set_y <- read.table(file =setY, header = TRUE,sep = separator,row.names = 1)

rownames(set_x) <- toupper(rownames(set_x))
rownames(set_y) <- toupper(rownames(set_y))

genes <- rownames(set_x)[rownames(set_x)%in% rownames(set_y)]


basename <- basename(setX)  
basename <- strsplit(split = ".csv", x = basename)[[1]][1]
path <- paste(dirname(setX),"Results", basename, sep = "/")  

clusteringNAME <- paste(path,"/",nclustX,"/",
                        basename,"_clustering.output.csv", sep = "")

write.table(x = set_x[genes,],paste(scrat_tmp.folder,"/",basename,"_X_Selected.csv",sep = ""),sep = separator, col.names = NA)

system(paste("cp ", clusteringNAME," ", paste(scrat_tmp.folder,"/",basename,"_X_clustering.output.csv",sep = ""),sep = ""))

new_clustering <- paste(basename,"_X_clustering.output.csv",sep = "")


Marker(clustering = new_clustering,
       path = markerX,
       pvalue = pvalueComet,
       scratch = scrat_tmp.folder, separator = separator)

### Ripeto procedimento per il set2
basename <- basename(setY)  
basename <- strsplit(split = "\\.", x = basename)[[1]][1]
path <- paste(dirname(setY),"Results", basename, sep = "/")  

clusteringNAME <- paste(path,"/",nclustY,"/",
                        basename,"_clustering.output.csv", sep = "")


write.table(x = set_y[genes,],paste(scrat_tmp.folder,"/",basename,"_Y_Selected.csv",sep = ""),
            sep = separator,
            col.names = NA)

new_clustering <- paste(basename,"_Y_clustering.output.csv",sep = "")

system(paste("cp ", clusteringNAME," ",scrat_tmp.folder,"/",new_clustering ,sep = ""))

clusteringNAME <- basename(clusteringNAME)

Marker(clustering = new_clustering,
       path = markerY,
       pvalue = pvalueComet,scratch = scrat_tmp.folder,
       separator = separator)


# system(paste("cp ",
#              X," ",
#              setY," ",
#              scrat_tmp.folder,"/",
#              sep="")
# )



# new_name_X
basename <- basename(setX)  
basename <- strsplit(split = "\\.", x = basename)[[1]][1]
nameX <- paste(scrat_tmp.folder,"/",basename,"_X_Selected.csv",sep = "")

# new_name_Y
basename <- basename(setY)  
basename <- strsplit(split = "\\.", x = basename)[[1]][1]
nameY <- paste(scrat_tmp.folder,"/",basename,"_Y_Selected.csv",sep = "")


variabili <- paste (basename(nameX), 
                    basename(nameY),
                    separator,
                    pvalueComet,
                    subdivision,
                    starting_genes,
                    contamination,
                    permutation,
                    threshold, 
                    sep = " ")



  #executing the docker job   VA CAMBIATO  CON IL COMANDO PER LANCIARE IL DOCKER
# params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,
#                 ":/scratch -d repbioinfo/bcsctemp Rscript /home/4.bis_New_BC_SCript.R ",
#                 basename(setX)," ", 
#                 basename(setAclustering)," ",
#                 basename(markerA), " ",
#                 basename(setB)," ", 
#                 basename(setBclustering)," ", 
#                 basename(markerB)," ",
#                 separator," ", 
#                 contamination," ",
#                 permutation," ",
#                 threshold, sep = "")
# 
    

params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,
                ":/scratch -d docker.io/repbioinfo/bcscswap.2021.01 Rscript /home/4.bis_New_BC_SCript.R ",
                variabili,
                sep = "")




############################ DA QUI NON TOCCO FINO A SEGNO ***
    
resultRun <- runDocker(group=group, params=params)

  #waiting for the end of the container work
  if(resultRun==0){
    #system(paste("cp ", scrat_tmp.folder, "/* ", data.folder, sep=""))
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
  system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/",
               substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))


  #Copy result folder
  cat("Copying Result Folder")
  system(paste("cp -r ",scrat_tmp.folder,"/BCsc ",data.folder,sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",
               data.folder, sep=""))
  setwd(home)
}

########################################################################################## ****

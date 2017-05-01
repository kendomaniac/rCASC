 
 
args <- commandArgs(TRUE)
if(length(args)==5){
system(paste("./2_analysis.sh",args[1]," ",args[2]," ",args[3]," ",args[4],args[5]))
}else{ system("./2_analysis.sh")}





 
 
args <- commandArgs(TRUE)
if(length(args)==3){
system(paste("./3_featureS.sh",args[1]," ",args[2]," ",args[3]))
}
if(length(args)==5){
system(paste("./3_featureS.sh",args[1]," ",args[2]," ",args[3]," ",args[4]," ",args[5]))
}
if(length(args)!=3 && length(args) != 5){
system("./3_featureS.sh")
}





 

 
args <- commandArgs(TRUE)
if(length(args)==7){
system(paste("./1_perm.sh",args[1]," ",args[2]," ",args[3]," ",args[4]," ",args[5]," ",args[6]," ",args[7] ))
}else{ system("./1_perm.sh")}





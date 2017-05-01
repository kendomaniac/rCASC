filename="path.txt"
while read -r line
do
    name="$line"
done < "$filename"
pathway=$name

if [ $# -ne 5 ] && [ $# -ne 3 ]
then


                    echo "Exists configuration file?y/n"
                    read answer2
                    if [ "$answer2" != "y" ] && [ "$answer2" != "n" ] 
                    then
                    echo "wrong letter"
                    exit
                    fi

                    if [ "$answer2" = "n" ] 
                    then 
                    echo "What is the matrix count name which you want analyze with feature selection?"
                    read matrixName
                    echo "MatrixName:"$matrixName > $pathway/options_3.txt

                    echo "What is the best number of cluster?"
                    read nCluster
                    echo "Number of Cluster:"$nCluster >> $pathway/options_3.txt

                    echo "Want you to use standard Z value?y/n"
                    read answer3
                                                    if [ "$answer3" != "y" ] && [ "$answer2" != "n" ] 
                                                    then
                                                    echo "wrong letter"
                                                    exit
                                                    exit
                                                    fi
                                                
                                                    if [ "$answer3" = "n" ] 
                                                    then 
                                                    echo "Set the initial value for Z sequence"
                                                    read firstZSeq
                                                    echo "firstZSeq:"$firstZSeq >> $pathway/options_3.txt
                                                    echo "Set the last value for Z sequence"
                                                    read lastZSeq
                                                    echo "lastZSeq:"$lastZSeq >> $pathway/options_3.txt
                                                
                                                    fi
                    echo "Want you to take genes specific for a D value or you want to collect all the characteristics genes in the D value range?answer y or n respectly"
                    read answer4 
                                                    if [ "$answer4" != "y" ] && [ "$answer4" != "n" ] 
                                                    then
                                                    echo "wrong letter"
                                                    exit
                                                    fi
                        echo "specificGenes?:"$answer4 >> $pathway/options_3.txt
                                                    
                    fi

fi

if [ $# -eq 5 ]
then
 echo "MatrixName:"$1 > $pathway/options_3.txt
echo "Number of Cluster:"$2 >> $pathway/options_3.txt
 echo "firstZSeq:"$3 >> $pathway/options_3.txt
 echo "lastZSeq:"$4 >> $pathway/options_3.txt
 echo "specificGenes?:"$5 >> $pathway/options_3.txt
fi

if [ $# -eq 3 ]
 then
 echo "MatrixName:"$1 > $pathway/options_3.txt
echo "Number of Cluster:"$2 >> $pathway/options_3.txt
 echo "specificGenes?:"$5 >> $pathway/options_3.txt

fi



echo "Want you to run docker with sudo?y/n"
read answer5 
if [ "$answer5" != "y" ] && [ "$answer5" != "n" ] 
then
echo "wrong letter"
exit
fi

if [ "$answer5" = "y" ] 
then 
sudo docker pull rcaloger/casc

nohup sudo docker run -i -v $pathway:/scratch rcaloger/casc Rscript /home/CASC/IndexOfDispersion.R  &

fi

if [ "$answer5" = "n" ] 
then 
docker pull rcaloger/casc
nohup docker run -i -v $pathway:/scratch rcaloger/casc Rscript /home/CASC/IndexOfDispersion.R  &

fi
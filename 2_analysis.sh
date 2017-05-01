filename="path.txt"
while read -r line
do
    name="$line"
done < "$filename"
pathway=$name
if [ $# -ne 5 ]
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
                    echo "What is the minimum relationship coefficent to define that two permutation clustered the same?"
                    read totIdentity 
                    echo "totIdentity:"$totIdentity > $pathway/options_2.txt

                    echo "What is the minimum Jaccard coefficent to define that two cluster are the same?"
                    read clusterIdentity 
                    echo "Cluster Identity:"$clusterIdentity >> $pathway/options_2.txt

                    echo "What is the min number of cluster you want analyze?"
                    read minNCluster
                    echo "min Number of Cluster:"$minNCluster >> $pathway/options_2.txt

                    echo "What is the max number of cluster you want analyze?"
                    read maxNCluster
                    echo "max Number of Cluster:"$maxNCluster >> $pathway/options_2.txt                    
                    
                    echo "What is the matrix count name you want analyze?"
                    read matrix 
                    echo "matrixName:"$matrix >> $pathway/options_2.txt

                    fi

              
fi
if [ $# -eq 5 ]
then 
 echo "totIdentity:"$1 > $pathway/options_2.txt
echo "Cluster Identity:"$2 >> $pathway/options_2.txt
echo "min Number of Clulster:"$3 >> $pathway/options_2.txt
  echo "max Number of Cluster:"$4 >> $pathway/options_2.txt
  echo "matrixName:"$5 >> $pathway/options_2.txt

fi 





      echo "Want you to run docker with sudo?y/n"
                    read answer3 
                    if [ "$answer3" != "y" ] && [ "$answer3" != "n" ] 
                    then
                    echo "wrong letter"
                    exit
                    fi
                    
                    
if [ "$answer3" = "y" ] 
then 
sudo docker pull rcaloger/casc

nohup sudo docker run -i -v $pathway:/scratch rcaloger/casc Rscript /home/CASC/Analysis.R  &
fi

if [ "$answer3" = "n" ] 
then 
 docker pull rcaloger/casc

nohup docker run -i -v $pathway:/scratch rcaloger/casc Rscript /home/CASC/Analysis.R  &

fi
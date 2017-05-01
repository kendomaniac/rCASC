echo "Are the pathway setted?y/n"
read answer1 


if [ "$answer1" != "y" ] && [ "$answer1" != "n" ] 
then
echo "wrong letter"
exit
fi

if [ "$answer1" = "n" ] 
then 
echo "Where you want your scratch folder be?"
read pathway 
pwd=$(pwd)
echo $pathway/scratch > $pwd/path.txt
mkdir $pathway/scratch 
echo "Where are your matrix count?"
read pathway2 
cp -rf $pathway2 $pathway/scratch 
fi

filename="path.txt"
while read -r line
do
    name="$line"
done < "$filename"
pathway=$name

if [ $# -ne 7 ]
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
                echo "What is the name of matrix file you want analyze?"
                read matrixName
                echo "matrixName:"$matrixName > $pathway/options.txt
                echo "How many permutation you want?"
                read nPerm
                echo "Number of Permutation:"$nPerm >> $pathway/options.txt
                echo "How many process in parallel you want?"
                read permAtTime
                echo "Permutation at Time:"$permAtTime >> $pathway/options.txt
                echo "How much core you want to reserve for SIMLR parallelization?"
                read cores
                echo "cores:"$cores >> $pathway/options.txt
                echo "What percentage you want to leave from totally cells for bootstrap?Write only number from 0 to 100"
                read percentage
                echo "Percentage cell to be remove:"$percentage >> $pathway/options.txt
                echo "Set the minimum number of Clusters"
                read minNumberOfCluster 
                echo "Min number of Cluster:"$minNumberOfCluster >> $pathway/options.txt
                echo "Set the max number of Clusters"
                read maxNumberOfCluster 
                echo "Max Number of Cluster:"$maxNumberOfCluster >> $pathway/options.txt
                fi

fi
if [ $# -eq 7 ]
then 
 echo "matrixName:"$1 > $pathway/options.txt
echo "Number of Permutation:"$2 >> $pathway/options.txt
echo "Permutation at Time:"$3 >> $pathway/options.txt
  echo "cores:"$4 >> $pathway/options.txt
echo "Percentage cell to be remove:"$5 >> $pathway/options.txt
 echo "Min number of Cluster:"$6 >> $pathway/options.txt
  echo "Max Number of Cluster:"$7 >> $pathway/options.txt
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
 nohup sudo docker run -i -v $pathway:/scratch rcaloger/casc Rscript /home/CASC/Main.R  &
fi

if [ "$answer3" = "n" ] 
then 
docker pull rcaloger/casc
nohup docker run -i -v $pathway:/scratch rcaloger/casc Rscript /home/CASC/Main.R  &

fi
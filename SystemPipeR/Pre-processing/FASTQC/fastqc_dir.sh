
if [ -z $1 ]
then
    	echo "sh fastqc_dir.sh directory output"
        exit -1
fi

dir=$1
shift


scriptsdir="/bigdata/messaoudilab/arivera/Scripts/fastqc"
directory=$(readlink -f "$dir")

for f in `$directory`*.gz
do
	echo "$f" 
 	out=$(dirname $f)
        sbatch --wrap "$scriptsdir/fastqc.sh $f "

done




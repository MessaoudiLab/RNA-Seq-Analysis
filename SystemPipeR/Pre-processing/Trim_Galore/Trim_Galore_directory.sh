#Takes a directory as input and loops through the files
if [ -z $1 ]
then
    	echo "sh trim_galore.sh directory score length"
        exit -1
fi

dir=$1
shift
qual=$1
shift
length=$1
shift

scriptsdir="/bigdata/messaoudilab/arivera/Scripts/trimgalore"
directory=$(readlink -f "$dir")

for f in `$directory`*.gz
do
  	out=$(dirname $f)
        echo "$f"
        sbatch --wrap "$scriptsdir/trim_galore_new.sh $f $qual $length "

done

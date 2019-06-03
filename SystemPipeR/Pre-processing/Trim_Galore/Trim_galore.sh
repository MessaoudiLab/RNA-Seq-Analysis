if [ -z $1 ]
then
    	echo "sh trim_galore.sh fastq quality length"
        exit -1
fi

fq1=$1
shift
qual=$1
shift
length=$1
shift
trim=$1
#shift
#out=$1

SEQFILE="$fq1"
if [ ! -z $1 ] || [[ $1 =~ '^[^-]' ]]
then
    	fq2=$1
	SEQFILE="$fq1 $fq2"
        shift
fi

echo $SEQFILE

module load fastqc/0.11.3
module load trim_galore/0.4.2

out=$(dirname $fq1)
#echo $out

trim_galore --clip_R1 13 --three_prime_clip_R1 2 -o $out -q $qual -length $length -fastqc $fq1
#fastqc -o $out -f fastq $SEQFILE

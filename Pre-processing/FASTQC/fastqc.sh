if [ -z $1 ]
then
    	echo "sh fastqc.sh fastq"
        exit -1
fi

fq1=$1
shift

SEQFILE="$fq1"
if [ ! -z $1 ] || [[ $1 =~ '^[^-]' ]]
then
    	fq2=$1
	SEQFILE="$fq1 $fq2"
        shift
fi

echo $SEQFILE

module load fastqc/0.11.3

out=$(dirname $fq1)
#echo $out

fastqc -o $out -f fastq --extract $fq1

#!/bin/bash
quantification=rsem
threads=8
while getopts r:c:tqn:S:X:-:h option
do
    case "${option}" in
	h) echo -e "./processRNASeq.sh [OPTIONS] [fastq files] \nOPTIONS:\n-r [ pe or se ]\n-c [ rsem or hisat ]\n-t [ only perform trimming ]\n-q [ perform quality control ]\n-n [ index name for hisat or index for rsem ]\n-S [ HISAT splice file ]\n--gtf [GTF file]\n--threads [ Number of Threads ]"; exit 0;;
	-)case "${OPTARG}" in
	      gtf)
		  gtf="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 )); 
		  ;;
	      threads)
		  threads="${!OPTIND}"; OPTIND=$(( $OPTIND + 1));
		  ;;
	  esac;;
	r) read_state=${OPTARG};
	   ;;
	c) quantification=${OPTARG}; 
	   ;;
	t) trim_only="yes";
	   ;;
	q) q_ctrl="yes";
	   ;;
	n) species=${OPTARG}; 
	   ;;
	S) splice=${OPTARG}; 
	   ;;
    esac
done
optnum=$(( $OPTIND - 1 ))
shift $optnum
if [ ! `echo $read_state` ];then
    auto_read=true
fi
#echo $auto_read $gtf $threads $read_state $quantification $trim_only $q_ctrl $species $splice
#exit 0

mkdir trim
mkdir raw
mkdir bam
mkdir QCData
echo -e "#!/bin/bash\n#PBS -N trimFastq\n#PBS -q copperhead \n#PBS -l procs=$threads" >> pbs.sh 
echo -e '#PBS -l walltime=40:00:00' >> pbs.sh
echo -e '#PBS -o /dev/null' >> pbs.sh
echo -e '#PBS -e /dev/null' >> pbs.sh
echo -e 'SHORT_JOBID=`echo $PBS_JOBID |cut -d. -f1`' >> pbs.sh
echo -e 'OUTFILE=$PBS_O_WORKDIR/$PBS_JOBNAME-$SHORT_JOBID.out' >> pbs.sh
echo -e 'ERRFILE=$PBS_O_WORKDIR/$PBS_JOBNAME-$SHORT_JOBID.err' >> pbs.sh
echo -e 'exec 1>$OUTFILE 2>$ERRFILE' >> pbs.sh
echo -e 'cd $PBS_O_WORKDIR' >> pbs.sh
echo -e 'module load samtools' >> pbs.sh
cat pbs.sh > temp.sh

num=`echo "$@" | wc -w`
while [ $num -gt 0 ]
do
    length=`less $1 |head -2 | grep -v ^@ | wc -c`; length=$(( $length - 1 ))
    name=`echo $1 | cut -f 1 -d"_"`
    name2=`echo $2 | cut -f1 -d"_"`
    lane=`echo $1 | tr "_" "\n" | grep L00`
    lane2=`echo $2 | tr "_" "\n" | grep L00`
    readnum=$(( `zcat $1 | wc -l` / 4 ))
    if [ "$quantification" = "hisat" ]; then
	name=$name\_$lane # For seperate processing of same samples on different lanes (RSEM works only on merged fasta files)
	name2=$name2\_$lane2
    fi
    if $auto_read; then
	if [ "$name" = "$name2" ]; then
            read_state="pe"
	else
            read_state="se"
	fi
    fi
    mkdir $name ; mkdir QCData/$name
    #Pair End Reads Processing
    if [ "$read_state" = "pe" ]; then
       	echo "trim_galore --length 20 --paired --phred33 --quality 25 --gzip -o $name --retain_unpaired $1 $2; mv $1 raw; mv $2 raw; mv $name trim" >> pbs.sh
	if [ "$trim_only" = "yes" ]; then
	    echo "exit 0" >> pbs.sh
	fi
	if [ "$quantification" = "rsem" ]; then
	    echo 'read1=`'"echo trim/$name/*val* | cut -f1 -d' '"'`' >> pbs.sh
	    echo 'read2=`'"echo trim/$name/*val* | cut -f2 -d' '"'`' >> pbs.sh
	    echo "rsem-calculate-expression --paired-end --star --star-gzipped-read-file --star-output-genome-bam --calc-ci --seed 9 --single-cell-prior -p $threads" '$read1 $read2' "$species $name" >> pbs.sh
	    echo "samtools sort -@ $threads $name.STAR.genome.bam -o $name\_sorted.bam; rm $name.STAR.genome.bam; mv $name\_sorted.bam ./bam" >> pbs.sh
	elif [ "$quantification" = "hisat" ]; then
	    echo "hisat2 -p $threads --dta-cufflinks -q --phred33 --known-splicesite-infile $splice -x $species -1 $read1 -2 $read2 -S $name.sam;" >> pbs.sh
	    echo "samtools view -@ $threads -b $name.sam -o $name.bam; rm $name.sam;" >> pbs.sh
	    echo "samtools sort -@ $threads $name.bam -o $name\_sorted.bam; rm $name.bam; mv $name\_sorted.bam ./bam" >> pbs.sh
       fi
       if [ "$q_ctrl" = "yes" ]; then
		echo "java -jar $QoRTs QC --minMAPQ 50 --maxReadLength $length --addFunctions calcDetailedGeneCounts --seqReadCt $readnum ./bam/$name\_sorted.bam $gtf ./QCData/$name" >> pbs.sh 
       fi
	qsub -N $name\_pbs pbs.sh
	tail -8 pbs.sh  
	cat temp.sh > pbs.sh
	shift 2
	num=$(($num-2))

    #Single End Read Processing
    elif [ "$read_state" = "se" ]; then
	    echo "trim_galore --length 20 --phred33 --quality 25 --gzip -o $name --fastqc  $1; mv $1 raw; mv $name trim;" >> pbs.sh
	    if [ "$quantification" = "rsem" ]; then
		echo 'read1=`'"echo ./trim/$name/*trimmed*fq*gz | cut -f1 -d' '"'`' >> pbs.sh
		echo "rsem-calculate-expression  --star --star-gzipped-read-file --calc-ci --fragment-length-mean 200 --fragment-length-sd 40 --seed 9 --single-cell-prior -p $threads" '$read1' "$species $name" >> pbs.sh
		echo "samtools sort -@ $threads $name.STAR.genome.bam -o $name\_sorted.bam; rm $name.STAR.genome.bam; mv $name\_sorted.bam ./bam" >> pbs.sh
	    elif [ "$quantification" = "hisat" ]; then
		name=$name\_$lane # For seperate processing of same samples on different lanes (RSEM works only on merged fasta files)
		echo "hisat2 -p $threads --dta-cufflinks -q --phred33 --known-splicesite-infile $splice -x $species -U $1 -S $name.sam;" >>pbs.sh
		echo "samtools view -@ $threads -b $name.sam -o $name.bam; rm $name.sam;" >> pbs.sh
		echo "samtools sort -@ $threads $name.bam -o $name\_sorted.bam; rm $name.bam; mv $name\_sorted.bam ./bam" >> pbs.sh
	    fi
	    if [ "$q_ctrl" = "yes" ]; then
		    echo "java -jar $QoRTs QC --minMAPQ 50 --maxReadLength $length --addFunctions calcDetailedGeneCounts --singleEnded --seqReadCt $readnum ./bam/$name\_sorted.bam $gtf ./QCData/$name" >> pbs.sh     
	    fi
	    
	qsub -N $name\_pbs pbs.sh
	tail -8 pbs.sh    
	cat temp.sh > pbs.sh;
	shift 1;
	num=$(($num-1))	
    fi
done
#mv *.bam bam
rm pbs.sh
rm temp.sh

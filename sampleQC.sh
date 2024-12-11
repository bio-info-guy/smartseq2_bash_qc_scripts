echo -e "#!/bin/bash\n#PBS -N trimFastq\n#PBS -q copperhead \n#PBS -l procs=$threads,mem=100GB" > pbs.sh
echo -e '#PBS -l walltime=5:00:00' >> pbs.sh
echo -e '#PBS -o /dev/null' >> pbs.sh
echo -e '#PBS -e /dev/null' >> pbs.sh
echo -e 'SHORT_JOBID=`echo $PBS_JOBID |cut -d. -f1`' >> pbs.sh
echo -e 'OUTFILE=$PBS_O_WORKDIR/$PBS_JOBNAME-$SHORT_JOBID.out' >> pbs.sh
echo -e 'ERRFILE=$PBS_O_WORKDIR/$PBS_JOBNAME-$SHORT_JOBID.err' >> pbs.sh
echo -e 'exec 1>$OUTFILE 2>$ERRFILE' >> pbs.sh
echo -e 'cd $PBS_O_WORKDIR' >> pbs.sh
echo -e 'module load samtools' >> pbs.sh
cat pbs.sh > temp.sh
gtf=`ls *gtf`
num=`echo "$@" | wc -w`
num=$(($num-1))
fname=$1; shift 1;
while [ $num -gt 0 ]
do
    name=`echo $1 | cut -f1 -d"_"`
    name2=`echo $2 | cut -f1 -d"_"`
    if [ "$name" = "$name2" ]; then
	echo -e "samtools merge -@8 -n $name\.bam $1 $2" >> pbs.sh
	echo -e "samtools sort -n -@8 -o $name.sam $name\.bam; rm $name.bam" >> pbs.sh
	pair=`samtools view -c -f 1 $1`
	num=$(($num-2))
	shift 2;
    else
	name_full=`echo $1 | cut -f1 -d "."`
	echo -e "samtools sort -n -@8 -o $name.sam $name_full\.bam;" >> pbs.sh
	pair=`samtools view -c -f 1 $1`
	num=$(($num-1))
	shift 1;
    fi
    if [ $pair != 0 ]; then
	echo -e "python sampleQC.py -n $name --sam $name.sam --gtf $gtf --rsem --rsem_ct $fname.rsem.ct --rsem_len $fname.eff.len --rsem_tpm $fname.tpm --sampling 1000 --spike " >> pbs.sh
	echo -e "python sampleQC.py -n $name --sam $name.sam --gtf $gtf --rsem --rsem_ct $fname.rsem.ct --rsem_len $fname.len --rsem_tpm $fname.tpm --sampling 1000000" >> pbs.sh
	echo -e "rm $name.sam;" >> pbs.sh
    else
	echo -e "python sampleQC.py -n $name --sam $name.sam --gtf $gtf --rsem --rsem_ct $fname.rsem.ct --rsem_len $fname.eff.len --rsem_tpm $fname.tpm --sampling 1000 --spike --se" >> pbs.sh
	echo -e "python sampleQC.py -n $name --sam $name.sam --gtf $gtf --rsem --rsem_ct $fname.rsem.ct --rsem_len $fname.eff.len --rsem_tpm $fname.tpm --sampling 1000000 --se" >> pbs.sh
	echo -e "rm $name.sam;" >> pbs.sh
    fi
    #cat pbs.sh;
    qsub -N $name\_customQC_pbs pbs.sh
    cat temp.sh > pbs.sh
done
rm temp.sh

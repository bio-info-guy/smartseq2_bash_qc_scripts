#!/usr/bin/bash
# use in QCData directory: bash mergeGeneExpressionQORTS.sh species *directories
name=$1; shift 1;

less $1/QC.geneCounts.detailed.txt.gz | cut -f1 | grep -v gene_id | sort -k 1b,1 > $name.qorts.ct
less $1/QC.geneCounts.detailed.txt.gz | cut -f1 | grep -v gene_id | sort -k 1b,1 > $name.qorts.cds.ct
less $1/QC.geneCounts.detailed.txt.gz | cut -f1 | grep -v gene_id | sort -k 1b,1 > $name.qorts.utr.ct
echo  >header.txt
for s in `echo $@`
do
#echo $s | rev | cut -f1 -d"/" | rev >> header.txt
echo $s >> header.txt
join -1 1 -2 1 <(cat $name.qorts.ct) <(less $s/QC.geneCounts.detailed.txt.gz | cut -f1,2 | grep -v GENEID | sort -k 1b,1) > $name.temp
cat $name.temp > $name.qorts.ct
join -1 1 -2 1 <(cat $name.qorts.cds.ct) <(less $s/QC.geneCounts.detailed.txt.gz | cut -f1,3 | grep -v GENEID  | sort -k 1b,1) > $name.temp
cat $name.temp > $name.qorts.cds.ct
join -1 1 -2 1 <(cat $name.qorts.utr.ct) <(less $s/QC.geneCounts.detailed.txt.gz | cut -f1,4 | grep -v GENEID  | sort -k 1b,1) > $name.temp
cat $name.temp > $name.qorts.utr.ct

done

cat header.txt | tr "\n" "\t" > $name.temp
echo  >> $name.temp
cat $name.qorts.utr.ct >> $name.temp
cat $name.temp |  tr " " "\t" > $name.qorts.utr.ct

cat header.txt | tr "\n" "\t"  > $name.temp
echo  >> $name.temp
cat $name.qorts.cds.ct >> $name.temp
cat $name.temp |  tr " " "\t" > $name.qorts.cds.ct

cat header.txt | tr "\n" "\t" > $name.temp
echo  >> $name.temp
cat $name.qorts.ct >> $name.temp
cat $name.temp | tr " " "\t" > $name.qorts.ct

rm header.txt
rm $name.temp

#!/usr/bin/bash
num=`echo "$@" | wc -w`
num=$(($num-1))
name=$1; shift 1;
cut -f1 $1 | grep -v gene_id > $name.eff.len
cut -f1 $1 | grep -v gene_id  > $name.rsem.ct
cut -f1 $1 | grep -v gene_id  > $name.tpm
echo  >header.txt
while [ $num -gt 0 ]
do
s_name=`echo $1 | cut -f1 -d"."`
echo $s_name >> header.txt
join -1 1 -2 1 <(cat $name.eff.len) <(cut -f1,4 $1| grep -v gene_id  ) > $name.temp
cat $name.temp > $name.eff.len
join -1 1 -2 1 <(cat $name.rsem.ct) <(cut -f1,5 $1| grep -v gene_id  ) > $name.temp
cat $name.temp > $name.rsem.ct
join -1 1 -2 1 <(cat $name.tpm) <(cut -f1,6 $1 | grep -v gene_id  ) > $name.temp
cat $name.temp > $name.tpm

num=$(($num-1))
shift 1;
done

cat header.txt | tr "\n" "\t" > $name.temp
echo  >> $name.temp
cat $name.eff.len >> $name.temp
cat $name.temp |  tr " " "\t" > $name.eff.len

cat header.txt | tr "\n" "\t"  > $name.temp
echo  >> $name.temp
cat $name.rsem.ct >> $name.temp
cat $name.temp |  tr " " "\t" > $name.rsem.ct

cat header.txt | tr "\n" "\t" > $name.temp
echo  >> $name.temp
cat $name.tpm >> $name.temp
cat $name.temp | tr " " "\t" > $name.tpm

rm header.txt
rm $name.temp

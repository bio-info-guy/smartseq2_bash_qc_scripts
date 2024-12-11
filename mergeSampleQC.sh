#!/usr/bin/bash
name=$1; shift 1;
echo -e  "sample_name\tnum_read\tspike_rate\tuniq_nospike\tsensitivity\tcomplexity\tevenness\tgap"> $name.stat
echo -e  "sample_name\tsp_read\tsp_uniq\tsensitivity\tsp_complexity\tsp_evenness\tsp_gap"> $name.spike.stat
#echo >> $name.stat; echo >> $name.spike.stat
for s in `echo $@`
do
    if [[ -f `echo $s/$s\_*spike*.stat` ]];then\
	cat $s/$s\_*spike*.stat >> $name.spike.stat
    else\
	echo -e "$s\t0\tNA\tNA\tNA\tNA\tNA" >> $name.spike.stat
    fi
    cat $s/$s\_*1000000*.stat >> $name.stat
done



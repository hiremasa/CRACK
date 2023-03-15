#!bin/bash

CT=../data/causal_tuebingen
CO=results/
count=1
while IFS=$'\t' read -r -a row
do
    echo "${row[0]} ${row[1]} ${row[6]}"
    if [ $count -eq 47 ] || [ $count -eq 70 ]
    then
        ./crack.run -s 1 -o ${CO}/tuebingen_ -x ${row[1]} -c -d ' ' -i ${CT}/${row[0]}.txt -a ${row[0]}
    else
        ./crack.run -s 1 -o ${CO}/tuebingen_ -x ${row[1]} -c -d ' ' -i ${CT}/${row[0]}.txt -a ${row[0]} -t 'i'
    fi
    count=$[$count + 1]
done < ${CT}/README_polished.tab

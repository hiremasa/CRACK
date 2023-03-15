#!bin/bash

CT=../data
CO=results/
while IFS=$'\t' read -r -a row
do
    echo "${row[0]} ${row[1]} ${row[5]} ${row[6]}"
    ./crack.run -s 3 -o ${CO}/MV_NCI_ -x ${row[1]} -d ' ' -c ${row[6]} -i ${CT}/${row[0]} -a ${row[0]}
done < ${CT}/mv_run.tab

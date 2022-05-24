#!/bin/bash

# This version works with abinit 9

qdir=$1

echo "000 TEMP" > totg_col.dat

# find volume
vol_line=$(grep -m 1 -r ucvol | head -1)
vol=$(echo $vol_line | awk '{print $(NF-1)}')

# assemble totg0's
echo "Assembling totg0's..."

for qfol in q*_at*; do 
    
    grep TOTG $qfol'/abinit.out' | awk '{$1="";print $0}' > common
    cat totg_col.dat common >> tmp
    rm common
    mv tmp totg_col.dat

done


# Sort
echo "Sorting..."

if (( qdir==1 )); then
    sort -gk1,1 -gk2,2 -gk3,3  totg_col.dat > tmp
    mv tmp totg_col.dat

elif (( qdir==2 )); then
    sort -gk1,1 -gk2,2 -gk4,4  totg_col.dat > tmp
    mv tmp totg_col.dat

elif (( qdir==3 )); then
    sort -gk1,1 -gk2,2 -gk5,5  totg_col.dat > tmp
    mv tmp totg_col.dat

else
    
    echo "No q direction specified, will sort from L to R"

    sort -gk1,1 -gk2,2 -gk3,3 -gk4,4 -gk5,5 totg_col.dat > tmp
    mv tmp totg_col.dat
fi

# include zion for calculation of Z
# Get ionic psp charges                                                                                                             
count=1
for ion in $(grep zion ./q1_at1/*.abo | awk '{print $3}'); do
    zion_type[$count]=$ion
    let count=count+1
done
count=1
for ion in $(grep " typat" ./q1_at1/*.abo | tail -n1); do
    if (( $count>1 )); then
        ct=$((count-1))
        zion[$ct]=${zion_type[$ion]}
    fi
    let count=count+1
done

natom=$(grep -w natom q1_at1/*.abo | tail -1 | awk '{print $2}')

# find number of lines
numlines=$(echo `wc -l totg_col.dat` | awk '{print $1}')
let numlines=numlines-1

replace_line=$(echo $numlines $vol $natom ${zion[@]})

sed -i "s/000 TEMP/$replace_line/g" totg_col.dat

# Sum over atoms
sum_dist_col.x

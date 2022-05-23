#!/bin/bash

qdir=$1

echo "000 TEMP" > totg_col.dat

# find volume
vol_line=$(grep -m 1 -r ucvol | head -1)
vol=$(echo $vol_line | awk '{print $(NF-1)}')

# assemble totg0's
echo "Assembling totg0's..."

for filein in $(find -name pindg.dat); do 

    sed -n "2,\$ p" $filein > common
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

# find number of lines
numlines=$(echo `wc -l totg_col.dat` | awk '{print $1}')
let numlines=numlines-1
sed -i "s/000 TEMP/$numlines $vol/g" totg_col.dat

# Sum over atoms
sum_dist_col.x

#!/bin/bash

file=oGaN.out

rprim_num=$(grep -m 1 -n rprim $file | awk '{print $1}')
rprim_num=${rprim_num%:*}

echo $rprim_num

count1=0
count2=0
count3=1
while [ count1 -lt 3 ]; do
    while [ count2 -lt 3 ]; do

	row=(( $rprim_num + count1 ))
	rprim[$count3]= $(sed "${NUM}q;d" $file | awk '{print $1}')

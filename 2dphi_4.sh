#!/bin/bash

# Extract dynamical matrix and BEC from output file of abinit calculation. Then run twoddyn_pi.x to take derivatives

# Cyrus Dreyer, Stony Brook and Flatiron CCQ


abo=$(ls *.abo)
if [ ! -z "$abo" ]; then
    fileout=$abo
else
    out=$(ls o*.out)
    if [ ! -z "$out" ]; then
	fileout=$out
    fi
fi

if [ -z "$fileout" ]; then
    echo 'ERROR: Cannot open output file'
    exit
fi

echo 'Output file:' $fileout

filein=*.in

# ***********************************
#    Read in lattice properties     *
#************************************

echo "2nd deriv mat dat" > dyn.dat

if [ ! -f bec.dat ]; then
    echo "Born effective charges" > bec.dat
    exist_bec=0
else
    bec_lines=$(wc -l bec.dat | awk '{print $1}')
    if (( bec_lines > 1 )); then
	echo "Will use existing bec.dat"
	exist_bec=1
    else
	echo "Born effective charges" > bec.dat
	exist_bec=0
    fi
fi


# size od dynamical matrix
nat=$(grep natom $filein | awk '{print $2}')
nblock=$(echo "$nat*3" | bc)
nqpt=$(grep -c "Dynamical matrix" $fileout)
echo "$nqpt $nat" > dyn.dat

# atom masses
ntypat=$(grep ntypat $filein | awk '{print $2}')
typat=$(grep -m 1 " typat" $fileout )
typat=${typat##*typat}
amu=$(grep -m 1 amu $fileout )
amu=${amu##*amu}
count1=0
for ii in $typat; do
    count=1
    for jj in $amu; do
	if [ $count -eq $ii ]; then
	    atmass[$count1]=$jj
	    let count1=count1+1
	    break
	fi
	let count=count+1
    done
done
echo >> dyn.dat
echo ${atmass[@]} >> dyn.dat
echo >> dyn.dat

# acel1
echo $(grep -m 1 acell $fileout | awk '{print $2" "$3" "$4}') >> dyn.dat

# ucvol
echo $(grep -m 1 ucvol $fileout | awk '{print $5}') >> dyn.dat

# rprim
line=$(grep -n rprim $filein)
nline=${line%%:*}
echo $(sed "${nline}q;d" $filein | awk '{print $2 " " $3 " " $4}') >> dyn.dat
let nline=nline+1
echo $(sed "${nline}q;d" $filein | awk '{print $1 " " $2 " " $3}') >> dyn.dat
let nline=nline+1
echo $(sed "${nline}q;d" $filein | awk '{print $1 " " $2 " " $3}') >> dyn.dat

echo >> dyn.dat

#xcart
# Not sure about this, for some reason no xcart in output for single atoms...
if  grep -q xcart $fileout; then 

    line=$(grep -m 1 -n xcart $fileout)
    nline=${line%%:*}

    count=0
    while [ $count -lt $nat ]; do
	
	if [ $count -eq 0 ]; then
	    echo $(sed "${nline}q;d" $fileout | awk '{print $2 " " $3 " " $4}') >> dyn.dat
	else
	    echo $(sed "${nline}q;d" $fileout | awk '{print $1 " " $2 " " $3}') >> dyn.dat
	fi
	
	let count=count+1
	let nline=nline+1
    done
else
    echo "0.0000000000 0.0000000000 0.0000000000" >> dyn.dat
fi

# get qpoints
count=1
count2=1
while [ $count2 -le $nqpt ]; do
    rfphon=$(grep -m $count rfphon $fileout |tail -1)
    rftst=$(echo $rfphon | awk '{print $2}')
    if [ $rftst -eq 1 ]; then
	rfphon=$(echo $rfphon | sed 's/rfphon//g')
	rfnum=$(echo $rfphon | awk '{print $1}')
	#qdat[$count2]=$rfnum
	qpt[$count2]=$(grep -m 1 " qpt$rfnum" $fileout | awk '{print $2 " " $3 " " $4}')
	let count2=count2+1
    fi
    let count=count+1
done

echo >> dyn.dat

#*********************************************
# Read in second deriv matrix from OUT file  *
#*********************************************

# Do this so we get it in cartesian coodrinates
count=1
countq=1
while [ $count -le $nqpt ]; do

    echo ${qpt[$count]}

    echo ${qpt[$count]} >> dyn.dat

    nline=$(grep -n -m $count "Dynamical matrix" $fileout | tail -n1)
    nline=${nline%%:*}
    let nline=nline+5

    count2=0
    while [ $count2 -lt $nblock ]; do
	count3=0
	while [ $count3 -lt $nblock ]; do
	    echo $(sed "${nline}q;d" $fileout) >> dyn.dat
	    let nline=nline+1
	    let count3=count3+1
	done
	let nline=nline+1
	let count2=count2+1
    done
    let count=count+1
    let countq=countq+1
done

#*********************************************
# Read in becs from OUT file  *
#*********************************************

if (( exist_bec==0 )); then
    nline=$(grep -n -m 2 "Effective charges" $fileout | tail -n1)
    nline=${nline%%:*}
    let nline=nline+6

    echo $nline

    count2=0
    while [ $count2 -lt $nblock ]; do
	count3=0
	while [ $count3 -lt 3 ]; do
	    echo $(sed "${nline}q;d" $fileout | awk '{print $5 " " $6}') >> bec.dat
	    let nline=nline+1
	    let count3=count3+1
	done
	let nline=nline+1
	let count2=count2+1
    done
fi

echo 'dyn.dat and bec.dat done'

# Get second derivatives with q	
twoddyn_pi.x


exit


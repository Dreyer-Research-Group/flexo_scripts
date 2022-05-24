#!/bin/bash

# Get the zero-displacement BECs
cd P0
dist_col_ab7.sh
cd q1_at1

# Get cartesian coordinates
nline=$(grep -n 'xcart' o*.out | tail -n1 | cut -d: -f1)
natom=$(grep natom o*.out | tail -n1 | awk '{print $2}')

count=1
while [ $count -le $natom ]; do

    line=$((nline+count-1))

    if (( $count==1 )); then
	cart_x[$count]=$(sed "${line}q;d" o*.out | awk '{print $2}')
	cart_y[$count]=$(sed "${line}q;d" o*.out | awk '{print $3}')
	cart_z[$count]=$(sed "${line}q;d" o*.out | awk '{print $4}')

    else
	cart_x[$count]=$(sed "${line}q;d" o*.out | awk '{print $1}')
	cart_y[$count]=$(sed "${line}q;d" o*.out | awk '{print $2}')
	cart_z[$count]=$(sed "${line}q;d" o*.out | awk '{print $3}')
    fi

    let count=count+1
done

# Get ionic psp charges
count=1
for ion in $(grep zion o*.out | awk '{print $3}'); do
    zion_type[$count]=$ion
    let count=count+1
done
count=1
for ion in $(grep " typat" o*.out | tail -n1); do
    if (( $count>1 )); then
	ct=$((count-1))
	zion[$ct]=${zion_type[$ion]}
    fi
    let count=count+1
done

cd ../

echo 'P^(0)_\alpha,\kappa\beta' > P0.dat
echo 'TEMP' >> P0.dat
echo 'kappa alpha beta         P^(0)                 cart_x              cart_y               cart_z (Bohr)' >> P0.dat
count=1
while read line; do
    if [ $count -gt 1 ]; then
	
	kappa=$(echo $line | awk '{print $1}')
	beta=$(echo $line | awk '{print $2}')
	p0[1]=$(echo $line | awk '{print $7}')
	p0[2]=$(echo $line | awk '{print $9}')
	p0[3]=$(echo $line | awk '{print $11}')

	idir=1
	while [ $idir -le 3 ]; do

	    # Add ionic part
	    if (( $idir==$beta )); then
		ze=$(echo ${p0[$idir]} | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
		bec=$(echo "scale=10; $ze+${zion[$kappa]}" | bc) 
	    else 
		bec=${p0[$idir]}
	    fi

	    echo '    '$kappa'     '$idir'    '$beta'     '$bec'     ' \
		${cart_x[$kappa]}'     '${cart_y[$kappa]}'     '${cart_z[$kappa]} >> P0.dat
	    let idir=idir+1
	done
    fi

    let count=count+1

done < totg_col.dat 

# Add number of lines to files for fortran input
P0_lines=$(wc -l P0.dat | awk '{print $1}' )
let P0_lines=P0_lines-3
sed -i 's/TEMP/'"$P0_lines"'/g' P0.dat

mv P0.dat ../


#cd $homedir # $homedir/


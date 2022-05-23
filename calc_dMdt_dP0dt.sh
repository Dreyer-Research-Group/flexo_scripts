#!/bin/bash

# This program will calculate the derivatives of M and P0
# right now \gamma=z

set -e

homedir=`pwd` # Run in Dz with Flexo.in in it

filein=Flexo.dat

fileout_dM=dMdt.dat
fileout_dZ=dZdt.dat

delt_r=$1

echo "Displacement:" $delt_r

echo "d M^z_k'b/dt_ka" > $fileout_dM
echo "TEMP" >> $fileout_dM
echo "kappa alpha kappa' gamma beta       dM/dt" >> $fileout_dM

echo "d Z^gamm_k'b/dt_ka" > $fileout_dZ
echo "TEMP" >> $fileout_dZ
echo "kappa alpha kappa' gamma beta       dZ/dt" >> $fileout_dZ


fxe_org=$(echo $homedir'/'$filein)

# Read in stuff from Flexo.dat
name=$(grep name $filein | awk '{print $2}')
nqpt=$(grep nqpt $filein | awk '{print $2}') 
atmin=$(grep rfatpol $filein | awk '{print $2}') 
atmax=$(grep rfatpol $filein | awk '{print $3}')  
nat=$((atmax-atmin+1))

# Get del_t. For now ONLY ORTHORHOMBIC!
a_times=$(grep acell $filein | awk  '{print $2}')
if [[ $a_times == *'*'* ]]; then
    acell[1]=$(echo $a_times | awk -F* '{print $2}')
    acell[2]=${acell[1]}
    acell[3]=${acell[1]}
else
    acell[1]=$(grep acell $filein | awk  '{print $2}')
    acell[2]=$(grep acell $filein | awk  '{print $3}')
    acell[3]=$(grep acell $filein | awk  '{print $4}')
fi

#TEST: get rid of acell??
#echo "TEST, NO ACELL!!!"
#acell[1]=1
#acell[2]=1
#acell[3]=1

#echo ${acell[@]}
#exit

del_t[1]=$(echo "scale=10; ${acell[1]}*$delt_r" | bc)
del_t[2]=$(echo "scale=10; ${acell[2]}*$delt_r" | bc)
del_t[3]=$(echo "scale=10; ${acell[3]}*$delt_r" | bc)


iat=1
while [ $iat -le $nat ]; do
    idir=1
    while [ $idir -le 3 ]; do

	folder=$(echo 'at_'$iat'_dir_'$idir)
	
	cd $folder # $homedir/at_i_dir_j
	echo "Entering "$folder
	
	# Read in Makb's and Z_akb
	cd "m_"$delt_r # $homedir/at_i_dir_j/m_0.001

	im=1
	while read m_in; do
	    if (( $im > 2 )); then
		imi=$((im-2))

		m_kappa[$imi]=$(echo $m_in | awk '{print $1}')
		m_alpha[$imi]=$(echo $m_in | awk '{print $2}')
		m_beta[$imi]=$(echo $m_in | awk '{print $3}')
		m_Makb[$imi]=$(echo $m_in | awk '{print $4}')		
	    fi

	    let im=im+1
	done < Makb.dat
	max_im_m=$((im-3))

	cd P1
	im=1
	while read m_in; do
	    if (( $im > 2 )); then
		imi=$((im-2))
		m_zkappa[$imi]=$(echo $m_in | awk '{print $1}')
		m_zalpha[$imi]=$(echo $m_in | awk '{print $2}')
		m_zbeta[$imi]=$(echo $m_in | awk '{print $3}')
		m_Zakb[$imi]=$(echo $m_in | awk '{print $4}')		
	    fi

	    let im=im+1
	done < Z_akb.dat
	max_im_mz=$((im-3))

	cd "../../p_"$delt_r # $homedir/Dz/at_i_dir_j/p_0.001
	im=1
	while read m_in; do
	    if (( $im > 2 )); then
		imi=$((im-2))
		p_kappa[$imi]=$(echo $m_in | awk '{print $1}')
		p_alpha[$imi]=$(echo $m_in | awk '{print $2}')
		p_beta[$imi]=$(echo $m_in | awk '{print $3}')
		p_Makb[$imi]=$(echo $m_in | awk '{print $4}')		
	    fi

	    let im=im+1
	done < Makb.dat
	max_im_p=$((im-3))

	cd P1 # $homedir/Dz/at_i_dir_j/p_0.001/P1
	im=1
	while read m_in; do
	    if (( $im > 2 )); then
		imi=$((im-2))
		p_zkappa[$imi]=$(echo $m_in | awk '{print $1}')
		p_zalpha[$imi]=$(echo $m_in | awk '{print $2}')
		p_zbeta[$imi]=$(echo $m_in | awk '{print $3}')
		p_Zakb[$imi]=$(echo $m_in | awk '{print $4}')		
	    fi

	    let im=im+1
	done < Z_akb.dat
	max_im_pz=$((im-3))

	# Check that all the lengths are consistent
	if (( $max_im_m != $max_im_p )) || (( $max_im_mz != $max_im_pz )); then
	    echo "The derivative files are not consistent!!"
	    exit
	fi
	
	echo "Makb and Zakb read"	
	cd ../../../ # $homedir/

	# ****************************
	# Calc d M^z_k^prime b/dt_ka *
	# ****************************

	# a = $idir, k = $iat, k^\prime = m/p_kappa, b = beta
	im=1
	while [ $im -le $max_im_m ]; do
	    
	    # Make sure M files are in the same order, etc.

	    if (( ${m_kappa[$im]} != ${p_kappa[$im]} )) || (( ${m_alpha[$im]} != ${p_alpha[$im]} )) || (( ${m_beta[$im]} != ${p_beta[$im]} )); then
		echo "Makb files are not consisitent!!"
		exit
	    fi
	    
	    # Deal with scientific notation in bc
	    pM=$(echo ${p_Makb[$im]} | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
	    mM=$(echo ${m_Makb[$im]} | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')	    

	    # Take the derivative
	    dMdt=$(echo "scale=20; ($pM-$mM)/(2*${del_t[$idir]})" | bc)
	    
	    
	    echo '  '$iat'     '$idir'     '${m_kappa[$im]}'      '${m_alpha[$im]}'     '${m_beta[$im]}'     '$dMdt >> $fileout_dM

	    let im=im+1
	done # loop through emements of M_kab

	echo "DONE calc M^z_k^prime b/dt_ka"

	# *******************************************
	# Calculate/ Collect P^(0)_alpha\kappa\beta * 
	# and derivatives for atomic displacements  *
	# *******************************************
	im=1
	while [ $im -le $max_im_mz ]; do
	    
	    # Make sure Z files are in the same order, etc.

	    if (( ${m_zkappa[$im]} != ${p_zkappa[$im]} )) || (( ${m_zalpha[$im]} != ${p_zalpha[$im]} )) || (( ${m_zbeta[$im]} != ${p_zbeta[$im]} )); then
		echo "Makb files are not consisitent!!"
		exit
	    fi
	    
	    # Deal with scientific notation in bc
	    pZ=$(echo ${p_Zakb[$im]} | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
	    mZ=$(echo ${m_Zakb[$im]} | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')	    

	    # Take the derivative
	    dZdt=$(echo "scale=20; ($pZ-$mZ)/(2*${del_t[$idir]})" | bc)
	    
	    echo '  '$iat'     '$idir'     '${m_zkappa[$im]}'      '${m_zalpha[$im]}'     '${m_zbeta[$im]}'     '$dZdt >> $fileout_dZ

	    let im=im+1
	done # loop through emements of Z_kab

	echo "DONE calc M^z_k^prime b/dt_ka"
	echo ""

	let idir=idir+1
    done # idir

    let iat=iat+1
done # iat


# Add number of lines to files for fortran input
M_lines=$(wc -l $fileout_dM | awk '{print $1}' )
Z_lines=$(wc -l $fileout_dZ | awk '{print $1}' )

let M_lines=M_lines-3
let Z_lines=Z_lines-3

sed -i 's/TEMP/'"$M_lines"'/g' $fileout_dM
sed -i 's/TEMP/'"$Z_lines"'/g' $fileout_dZ

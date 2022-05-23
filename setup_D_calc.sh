#!/bin/bash

# This program will run the scripts necessary to get D^z after the
# calculations setup by setup_M_calc.sh are run
# Should be run in the Dz folder

set -e

delt_r=$1

homedir=`pwd`

for fol in at_*; do

    echo "entering "$fol

    cd $fol # $homedir/at_i_dir_j
    cd "m_"$delt_r # $homedir/at_i_dir_j/m_0.001
    
    count=1
    while [ $count -le 2 ]; do
	# get derivatives_piezo.dat
	cd P1 # $homedir/at_i_dir_j/m_0.001/P1
	echo "Running dist_col and flexo_derivatives in " `pwd`
	dist_col_2.sh     
	#flexo_derivatives.x
	onlyeBohr_flexo_deriv.x

	cd ../P2 # $homedir/at_i_dir_j/m_0.001/P2
	echo "Running dist_col and flexo_derivatives in " `pwd`
	dist_col_2.sh
	echo "Gathered totg0s"
	#flexo_derivatives.x
	onlyeBohr_flexo_deriv.x

	# get M_kab.dat files
	cd ../ #$homedir/at_i_dir_j/m_0.001/
	echo "Running calc_Makb in " `pwd`
	calc_Makb.x
	
	cd "../p_"$delt_r #$homedir/at_i_dir_j/p_0.001/
	let count=count+1
    done
    
    cd $homedir # $homedir/
done

echo "DONE setting up Makb and Zakb"

# Calculate the displacemnt derivatives
echo "Running calc_dMdt_dP0dt"
calc_dMdt_dP0dt.sh $delt_r
cd $homedir # $homedir/

# Get the zero-displacement BECs
echo "Getting the zero-displacement BECs"
get_P0.sh

# Finally calculate D tensor
echo "Calculate D tenosr"
calc_Dz_kakb.x

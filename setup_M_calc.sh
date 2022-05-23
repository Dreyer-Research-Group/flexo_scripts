#!/bin/bash

# This program will setup the calculation of the antisymetric dynamic quadrapoles

#set -e

# Enter the displacement
disp=$1
mdisp=$(echo "(-1)*$disp" | bc)


if [ -z "$disp" ]; then
    echo "No displacement specified"
    exit
fi

filein=Flexo.dat
homedir=`pwd`

fxe_org=$(echo $homedir'/'$filein)

# Read in stuff
name=$(grep name $filein | awk '{print $2}')
nqpt=$(grep nqpt $filein | awk '{print $2}') 
atmin=$(grep rfatpol $filein | awk '{print $2}') 
atmax=$(grep rfatpol $filein | awk '{print $3}')  
nat=$atmax-$atmin+1


# get magnitude of q points
delq=0.01
count=1
while [ $count -le $nqpt ]; do
    qpts[$count]=$(echo "scale=5; ($count-1)*$delq" | bc | awk '{printf "%f", $0}')
    let count=count+1
done

# Get reduced coord
natom=$(grep natom $filein | awk '{print $2}')
nxred=$(awk '/xred/{ print NR; exit }' $filein)

count=1
while [ $count -le $natom ]; do
    line=$((count+nxred))
    xredx[$count]=$(sed "${line}q;d" $filein | awk '{print $1}')
    xredy[$count]=$(sed "${line}q;d" $filein | awk '{print $2}')
    xredz[$count]=$(sed "${line}q;d" $filein | awk '{print $3}')

    let count=count+1
done


mkdir Dz
cp $filein Dz
cd Dz # $homedir/Dz

#**************************************
# Make directory for undisplaced BECs *
#**************************************
mkdir P0
cd P0 # $homedir/Dz/P0

# Copy other files here
cp $(echo $homedir'/*') .

# Just need to do first qpt, which should be gamma 
nqpt_line=$(grep nqpt $filein)

sed -i 's/nqpt 3/nqpt 1/g' $filein
sed -i 's/nqpt 9/nqpt 1/g' $filein

# Run flexo_abin_gen.sh     
flexo_abin_gen.sh

echo "Done with P0"

cd ../ # $homedir/Dz 

#*************************************
# Make directories for displacements *
#*************************************
count=1
while [ $count -le $natom ]; do
    dircount=1
    while [ $dircount -le 3 ]; do

	dir_name=$(echo 'at_'$count'_dir_'$dircount)
	mkdir $dir_name
	cd $dir_name # $homedir/Dz/at_x_dir_y
	
	# Plus and minus deformation
	pm_count=1
	while [ $pm_count -le 2 ]; do

	    # plus deformation
	    if (( $pm_count==1 )); then
		disp_pm=$disp
		mkdir "p_"$disp
		cd "p_"$disp # $homedir/Dz/at_x_dir_y/p_0.01
		
	    elif (( $pm_count==2 )); then
		disp_pm=$mdisp
		mkdir "m_"$disp
                cd "m_"$disp # $homedir/Dz/at_x_dir_y/m_0.01
	    fi

	    # directions of q vector (now x and y)
	    qdir=1
	    while [ $qdir -le 2 ]; do
		
		dir_name=$(echo 'P'$qdir)
		mkdir $dir_name
      		cd $dir_name # $homedir/Dz/at_x_dir_y/m_0.01/P
		
		# Copy other files here
		cp $(echo $homedir'/*.files') .

		
		fhi_test=$homedir'/*.fhi'
		psp8_test=$homedir'/*.psp8'
		if [ -f $(echo $fhi_test | awk '{print $1}') ]; then
		    cp $(echo $homedir'/*.fhi') .
		elif [ -f $(echo $psp8_test | awk '{print $1}') ]; then
	    	    cp $(echo $homedir'/*.psp8') .
		fi

		# Head of the Flexo.dat file:
		echo '# Flexo file for '$dir_name > Flexo.dat
		echo '' >> Flexo.dat
		echo 'name '$name >> Flexo.dat
		echo '' >> Flexo.dat
		echo 'nqpt '$nqpt >> Flexo.dat
		echo '' >> Flexo.dat
		echo 'rfatpol 1 '$natom >> Flexo.dat
		echo '' >> Flexo.dat
		echo 'qpts' >> Flexo.dat
		
		# Print q points to Flexo.dat
		qcount=1
		while [ $qcount -le $nqpt ]; do
		    
		    if (( $qdir == 1 )); then
			echo ${qpts[$qcount]} '0.0 0.0' >> Flexo.dat
		    elif (( $qdir == 2 )); then
			echo '0.0 ' ${qpts[$qcount]} ' 0.0' >> Flexo.dat
		    elif (( $qdir == 3 )); then
			echo '0.0  0.0 ' ${qpts[$qcount]} >> Flexo.dat
		    fi
		    
		    let qcount=qcount+1
		done

		# Lets get all of the directions done
		echo '' >> Flexo.dat
		echo 'rfdir 1 1 1' >> Flexo.dat
		
		# Right now specify the split
		#echo 'split '$nqpt' '$natom >> Flexo.dat
		echo 'split 1 1' >> Flexo.dat
		echo '' >> Flexo.dat

		# Now for the common input variables
		ncom=$(grep -n "Common input" $fxe_org)
		ncom=${ncom%%:*} 

		sed -n "$ncom,\$ p" $fxe_org > common	
		cat Flexo.dat common >> tmp
		mv tmp Flexo.dat
		rm common

		# Finally, add in atomic displacements
		nxred=$(awk '/xred/{ print NR; exit }' Flexo.dat)
		xred_mod_ln=$((nxred+count))
	   
		if (( $dircount==1 )); then
		    new_line_x=$(echo "scale=10; ${xredx[$count]}+($disp_pm)" | bc | awk '{printf "%f", $0}')
		    new_line=$(echo $new_line_x' '${xredy[$count]}' '${xredz[$count]})
		    sed -i "$xred_mod_ln"'s/.*/'"$new_line"'/g' Flexo.dat

		elif (( $dircount==2 )); then
		    new_line_y=$(echo "scale=10; ${xredy[$count]}+($disp_pm)" | bc | awk '{printf "%f", $0}')
		    new_line=$(echo ${xredx[$count]}' '$new_line_y' '${xredz[$count]})
		    sed -i "$xred_mod_ln"'s/.*/'"$new_line"'/g' Flexo.dat

		elif (( $dircount==3 )); then
		    new_line_z=$(echo "scale=10; ${xredz[$count]}+($disp_pm)" | bc | awk '{printf "%f", $0}')
		    new_line=$(echo ${xredx[$count]}' '${xredy[$count]}' '$new_line_z)
		    sed -i "$xred_mod_ln"'s/.*/'"$new_line"'/g' Flexo.dat
		fi

		# Run flexo_abin_gen.sh
		flexo_abin_gen.sh	       

		let qdir=qdir+1
		cd ../ # $homedir/Dz/at_x_dir_y/pm_0.01
	    done # qdir for direction of q vec
	    
	    cd ../ # $homedir/Dz/at_x_dir_y/ 
	    let pm_count=pm_count+1
	done # pm for plus and minus displacements
	
	let dircount=dircount+1
	cd ../ # $homedir/Dz/
    done # dircount for direction of atomic displacement

    let count=count+1
done # count for number of atoms



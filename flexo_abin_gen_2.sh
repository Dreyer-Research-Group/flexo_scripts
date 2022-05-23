#!/bin/bash

# This program will read from a Flexo.dat file and generate abinit input file(s) for calulation of CI flexo constants

filein=Flexo.dat

# Read in stuff
name=$(grep name $filein | awk '{print $2}')
nqpt=$(grep nqpt $filein | awk '{print $2}') 
atmin=$(grep rfatpol $filein | awk '{print $2}') 
atmax=$(grep rfatpol $filein | awk '{print $3}')  
nat=$atmax-$atmin+1

# Test if metric perturbation
if ( grep -q "metric" Flexo.dat ); then
    metric=1
else
    metric=0
fi

# Test if diamag suscept
if ( grep -q "diamag" Flexo.dat ); then
    diamag=1
else
    diamag=0
fi


rfdir=$(grep rfdir $filein | awk '{print $2 " " $3 " " $4}')
#rfdiry=$(grep rfdir $filein | awk '{print $3}')
#rfdirz=$(grep rfdir $filein | awk '{print $4}')

splitq=$(grep split $filein | awk '{print $2}')
splitat=$(grep split $filein | awk '{print $3}')
splitdir=$(grep split $filein | awk '{print $3}')

# Check directions

qline=$(grep -n qpts $filein)
qline=${qline%%:*} 

count=1
nqzero=0
while [ $count -le $nqpt ]; do
    let qline=qline+1
    qpt[$count]=$(sed "${qline}q;d" $filein)

    # test for q = gamma. SHOULD BE FIRST IF PRESENT!!!!!
    qx=$(echo ${qpt[$count]} | awk '{print $1}')
    qy=$(echo ${qpt[$count]} | awk '{print $2}') 
    qz=$(echo ${qpt[$count]} | awk '{print $3}') 

    qzero=$(echo "scale=10; ($qx^2+$qy^2+$qz^2) < 1*10^(-9)" | bc)
    
    if [[ $qzero -eq 1 ]]; then
	nqzero=$count
    fi
    let count=count+1
done


# how many loops?
qgrp=$(( ($nqpt+$splitq-1)/$splitq ))
atgrp=$(( ($nat+$splitat-1)/$splitat ))

# loop over splits
countq=1
qstart=1
while [ $countq -le $splitq ]; do
    countat=1
    atstart=$atmin
    while [ $countat -le $splitat ]; do

	# Make sure we're not done
	if [[ $atstart -gt $atmax ]]; then
	    break
	fi

	folname="q"$countq"_at"$countat
	mkdir $folname
	
	cp * $folname
	cd $folname
		
	# Get last element in group
	qfin=$(( $qstart+$qgrp-1 ))
	if [[ $qfin -gt $nqpt ]]; then
	    qfin=$nqpt
	fi
	qnum=$(echo "scale=0; $qfin-$qstart+1" | bc)

	atfin=$(( $atstart+$atgrp-1 ))
	if [[ $atfin -gt $atmax ]]; then
	    atfin=$atmax
	fi
	rfatpolgrp=$atstart" "$atfin

	# Number of datasets (less for diamag):
	if [ $diamag -eq 0 ]; then 
	    ndtset=$(echo "scale=0; 1+($qfin-($qstart)+1)*3" | bc)
	else
	    ndtset=$(echo "scale=0; 1+($qfin-($qstart)+1)*2" | bc)
	fi

	qzero=0
	if [[ $nqzero -le $qfin ]] && [[ $nqzero -ge $qstart ]]; then
	    qzero=1
	    let ndtset=ndtset-1
	fi

	#******************
	# Make input file *
	#******************
	echo "# Flexo autogen input file" > $name
	echo " " >> $name
	echo "ndtset  $ndtset" >> $name
	echo " " >> $name

	# SCF Ground state calculation
	echo '#Ground state calculation ##############################################' >> $name
	echo 'kptopt1   1' >> $name
	echo 'tolvrs1   1.0d-18' >> $name
	echo 'prtden1   1   ' >> $name
	echo 'nqpt1   0' >> $name
	echo 'getwfk1   0 ' >> $name
	echo '########################################################################' >> $name
	echo " " >> $name


	# NSCF GS calculation
	echo '#Non-self consistent ground-state calculation##########################' >> $name
	echo 'getwfk   1' >> $name
	echo 'nqpt   1' >> $name
	echo 'symfxe 1' >> $name # No symmetry in looppert to be safe
	echo 'prepgkk 1' >> $name # Do all of the perturbations
	echo " " >> $name

	count=$qstart
	setnum=2
	while [ $count -le $qfin ]; do
	    
	    if [[ $nqzero -ne $count ]]; then
		echo "qpt"$setnum"  "${qpt[$count]} >> $name
		echo "getden"$setnum"  1" >> $name
		echo "tolwfr"$setnum"  1.0d-18" >> $name
		echo "iscf"$setnum"  -2" >> $name
		echo " " >> $name

		let setnum=setnum+1
	    fi
	    let count=count+1
	done

	# Static phonon calc
	
	if [ $diamag -eq 0 ]; then # Don't need phonon/metric for diamag

	    echo '#Static Response Function calculation####################################' >> $name
	    if [ $metric -eq 0 ]; then
		echo "rfatpol  "$rfatpolgrp >> $name
		echo "rfdir  "$rfdir >> $name
		echo " " >> $name

		count=$qstart
		while [ $count -le $qfin ]; do
		    echo "qpt"$setnum"  "${qpt[$count]} >> $name
		    echo "rfphon"$setnum"  1" >> $name

		    if [[ $nqzero -ne $count ]]; then
			wfq=$(( $setnum-$qnum ))
			echo "getwfq"$setnum"  "$wfq >> $name
		    fi

		    echo "nogzero"$setnum"  1" >> $name
		    echo "tolvrs"$setnum"  1.0d-8" >> $name
		    echo " " >> $name
		    let count=count+1
		    let setnum=setnum+1
		done
	    else
		# Static metric calc
		echo "rfatpol  "$rfatpolgrp >> $name
		echo "rfdir  "$rfdir >> $name
		echo " " >> $name

		count=$qstart
		while [ $count -le $qfin ]; do
		    echo "qpt"$setnum"  "${qpt[$count]} >> $name
		    echo "rfuser"$setnum"  2" >> $name

		    if [[ $nqzero -ne $count ]]; then
			wfq=$(( $setnum-$qnum ))
			echo "getwfq"$setnum"  "$wfq >> $name
		    fi

		    echo "nogzero"$setnum"  1" >> $name
		    echo "useylm"$setnum"  1" >> $name
		    echo "tolvrs"$setnum"  1.0d-8" >> $name
		    echo " " >> $name
		    let count=count+1
		    let setnum=setnum+1
		done
	    fi # metric
	fi # diamag

	# Adiabatic response calc
	echo '#Adiabatic Resonse Function############################################' >> $name
	
	# For diamag, had skipped this above
	if [ $diamag -eq 1 ]; then
	    echo "rfatpol  "$rfatpolgrp >> $name
	    echo "rfdir  "$rfdir >> $name
	    echo " " >> $name
	fi

	count=$qstart
	while [ $count -le $qfin ]; do
	    echo "qpt"$setnum"  "${qpt[$count]} >> $name
	    echo "rfphon"$setnum"  1" >> $name
	    echo "iscf"$setnum"  -3" >> $name

	    # Don't need this for diamag
	    if [ $diamag -eq 0 ]; then
		wfone=$(( $setnum-$qnum ))
		echo "get1wf"$setnum"  "$wfone >> $name
	    fi

	    if [[ $nqzero -ne $count ]]; then
		
		# Slightly different for diamag
		if [ $diamag -ne 1 ]; then
		    wfq=$(( $setnum-2*$qnum ))
		else
		    wfq=$(( $setnum-$qnum ))
		fi
		echo "getwfq"$setnum"  "$wfq >> $name
	    fi
	    
	    # Trigger diamag or adiabatic calc
	    if [ $diamag -eq 1 ]; then
		echo "adcalc"$setnum"  2" >> $name
	    else
		echo "adcalc"$setnum"  1" >> $name
	    fi

	    echo "tolwfr"$setnum"  1.0d-16" >> $name
	    
	    if [ $metric -eq 1 ]; then 
		echo "metcalc"$setnum"  1" >> $name
	    fi

	    echo " " >> $name
	    let count=count+1
	    let setnum=setnum+1
	done

	atstart=$(( $atfin+1 ))
	let countat=countat+1

	# Add the Common input variables
	ncom=$(grep -n "Common input" ../$filein)
	ncom=${ncom%%:*} 

	sed -n "$ncom,\$ p" ../$filein > common	
	cat $name common >> tmp
	mv tmp $name
	rm common

	cd ..

    done

    qstart=$(( $qfin+1 ))
    let countq=countq+1
done


# Create folders for split


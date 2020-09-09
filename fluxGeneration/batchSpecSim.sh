#!/bin/bash

INPUT=$CRPROPA_DIR/CRPropa3/build # Your CRPropa dir here
cd $INPUT
#EXPONENTS="20.1 21.1 23.1" # Emax
EXPONENTS="20.0 21.0 23.0" # Emax
SEZS="1.5"
#SEZS="1.5"
ALPHAS="2.0"
CS="1"
SS="1"

NUMBER="10000000" # number of protons to propagate from the source
BATCH="1" # always on for this script

for sourceEvolutionRedshiftCap in `echo $SEZS`; do
    for sourceEvo in `echo $SS`; do
	for cutoff in `echo $CS`; do
	    for alpha in `echo $ALPHAS`; do
		for exp in `echo $EXPONENTS`; do
		    fileName=$INPUT/cosmicRayMaxEnergyComparison"_E"$exp"_A"$alpha"_c"$cutoff"_s"$sourceEvo"_z"$sourceEvolutionRedshiftCap.sh
		    cd $(pwd) > $fileName
		    #echo source $HOME/anitaBuildTool/components/icemc/env.sh >> $fileName
		    echo source $HOME/.bashrc >> $fileName
		    echo cd $INPUT >> $fileName
		    echo source $CRPROPA_DIR\"/bin/activate\" >> $fileName
		    echo source /opt/rh/devtoolset-7/enable >> $fileName
		    echo source /opt/rh/python27/enable >> $fileName
		    echo source /cvmfs/sft.cern.ch/lcg/contrib/gcc/8.2.0/x86_64-centos7-gcc8-opt/setup.sh >> $fileName
		    echo python $INPUT/simulatedNeutrinoFlux.py -e $exp -a $alpha -c $cutoff -s $sourceEvo -z $sourceEvolutionRedshiftCap -n $NUMBER -b $BATCH >> $fileName
		    chmod +x $fileName
		    echo "*****************************************************************************"
		    echo "Submitting script $fileName to batch farm..."
		    pwd
		    cat $fileName
		    echo "*****************************************************************************"
		    #qsub -l mem=16gb:nodes=1:ppn=8:typea -q longc7 $fileName
		    #qsub -l mem=16gb:nodes=1:ppn=8:typea -q mediumc7 $fileName
		    #qsub -l mem=16gb:nodes=1:ppn=8:typea -q mediumc7 $fileName		THESE ARE BATCH
		    qsub -l mem=16gb:nodes=1:ppn=8:typea -q shortc7 $fileName		
		    #python $INPUT/simulatedNeutrinoFlux.py -e $exp -a $alpha -c $cutoff -s $sourceEvo -z $sourceEvolutionRedshiftCap -n $NUMBER -b $BATCH  # NON BATCH, USE WHEN BATCH FARM IS FULL
		    # request more mem, individual node, multithreaded, long q
		    echo "*****************************************************************************"
		    echo ''
		done
	    done
	done
    done
done

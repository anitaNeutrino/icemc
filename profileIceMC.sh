## Creates an output profiler for iceMC for a given number of events, and submits to a batch farm
## Input usage: ./profileICEMC.sh <NUMBER_OF_EVENTS>
## Output usage: kcachegrind iceMC_profile_out.<NUMBER_OF_EVENTS>
#!/bin/bash

## Customise this for where you wish to look
INPUT=/home/batten/anitaBuildTool/components/icemc
#OUTPUT=/unix/anita4/berg/testOutput/profiling
OUTPUT=$INPUT/profiling

## IceMC vars
NUMEVENTS=$1
EXPONENT=20
RUN=7357

## Valgrind input
fileName=$INPUT/iceMC_profile_in.$NUMEVENTS.sh
cd $(pwd) > $fileName
echo source $INPUT/env.sh >> $fileName
echo cd $INPUT >> $fileName
echo source $INPUT/setup.sh >> $fileName
outdir=$OUTPUT/$NUMEVENTS
rm -rf $outdir
mkdir $outdir
## Valgrind output: a callgrind file to use with kcachegrind
echo valgrind --tool=callgrind --callgrind-out-file=$outdir/iceMC_profile_out.$NUMEVENTS ./icemc -i $INPUT/inputs.anita4.conf -n $NUMEVENTS -e $EXPONENT -r $RUN >> $fileName
chmod +x $fileName
echo "*****************************************************************************"
echo "Submitting script $fileName to batch farm..."
cat $fileName
echo "*****************************************************************************"
## Queue job
qsub -q longc7 -o ${outdir}/output.log -e ${outdir}/error.log $fileName 
echo "*****************************************************************************"
echo ''

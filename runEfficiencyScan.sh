#!/bin/bash

source ~/ANITA/anita3/anitaBuildTool/env.sh
source ~/ANITA/anita3/anitaBuildTool/components/icemc/setup.sh

inputFile=inputs.anita3_trigEffScan.conf

for att in `seq -55 -40`;do

    sed -i.bak '/Off-axis attenuation/d' ${inputFile}
    echo "Off-axis attenuation: 0, 0, $att, $att, 0 # Attenuation applied to central phi sectors and the two adjecent ones (when 999 no signal in those antennas)" >> ${inputFile}

    att=$(($att * -1))
    output=efficiencyScanPulseAtAntennaNewTriggerBenRotter-7dBNewNoise/Att${att}db
    mkdir -p ${output}
    echo $output
    ./testInputAfterAntenna -i ${inputFile} -o ${output} -r 1 -n 100    

done

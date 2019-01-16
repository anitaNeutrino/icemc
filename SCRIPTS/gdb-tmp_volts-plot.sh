#!/usr/bin/env bash

GAIN_VOLTS_COL_FILE=$1
TRIG_VOLTS_COL_FILE=$2
GAIN_VS_TRIGGER_VOLTS=$3

TMP_GNUPLOT_SCRIPT=GDB/amp_eq_1/_gdb-tmp_volts-plot.gpl
cat << EOF  > $TMP_GNUPLOT_SCRIPT
set out "$GAIN_VS_TRIGGER_VOLTS"
set term post color solid 18
# plot [200: 300]
plot [:]\
       	"$GAIN_VOLTS_COL_FILE" with lines lt 1 lw 3,\
       	"$TRIG_VOLTS_COL_FILE" with lines lt 2
EOF

gnuplot --persist $TMP_GNUPLOT_SCRIPT

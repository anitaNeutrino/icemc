#!/usr/bin/env bash

GAIN_VOLTS_COL_FILE=$1
GAIN_VS_TRIGGER_VOLTS=$2

TMP_GNUPLOT_SCRIPT=GDB/amp_eq_1/_gdb-tmp_volts-plot.gpl
cat << EOF  > $TMP_GNUPLOT_SCRIPT
set out "$GAIN_VS_TRIGGER_VOLTS"
set term post color solid 18
plot "$GAIN_VOLTS_COL_FILE" with lines
EOF

gnuplot --persist $TMP_GNUPLOT_SCRIPT

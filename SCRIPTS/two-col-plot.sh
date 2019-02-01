#!/usr/bin/env bash

FILE1=$1
FILE2=$2
OUT=$3
XRANGE=$4
YRANGE=$5
SCALE=$6
SHIFT=$7
TMP_GNUPLOT_SCRIPT=$8

echo "Temporary gnuplot script is: $TMP_GNUPLOT_SCRIPT"

cat << EOF  > "$TMP_GNUPLOT_SCRIPT"
set out "$OUT"
set term post color solid 18
plot $XRANGE $YRANGE \
       	"$FILE1" using 0:(column(1) * 1.0) with lines  lt 1 lw 3,\
       	"$FILE2" using (column(0) $SHIFT):(column(1) * $SCALE) with lines lt 2
EOF

gnuplot $TMP_GNUPLOT_SCRIPT

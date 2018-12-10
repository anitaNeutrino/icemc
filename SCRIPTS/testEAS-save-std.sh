#!/usr/bin/env bash

EXE=$1
INP1=$2
INP2=$3
OUT=$4

$EXE -i $INP1 -o bvv_test -r 123 -n 1 -s $INP2 2>&1 | tee $OUT

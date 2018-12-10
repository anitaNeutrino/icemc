#!/usr/bin/env bash

INP=$1
PREFIX=$2
OUT=$3
TAB=$(echo $'\t')

grep "$PREFIX" $INP | sed -e "s/$PREFIX/$TAB/g" > $OUT

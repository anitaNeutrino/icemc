#!/usr/bin/env bash
INPUT_EPS_FILE="$1"

ghostscript -c "<</Orientation 3>> setpagedevice" -f $INPUT_EPS_FILE -c quit



#!/usr/bin/env bash

gawk '
  BEGIN{data=0}
  /BIZDOOMCAR/ {data=1; next}
  { if (data) print }
' $1 | 
sed 's/.*{\(.*\)}$/\1/;s/, */\n/g' > $2

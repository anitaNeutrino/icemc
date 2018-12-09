#!/usr/bin/env bash

sed 's/.*{\(.*\)}$/\1/;s/, */\n/g' $1 > $2

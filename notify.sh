#!/usr/bin/env bash

# $1: executable name.
# $2: variable name.
# $3: new value of the variable.
pidof $1 || exit 1

gdb -n -q -batch -ex "attach $(pidof $1)" \
  -ex "set $2 = \"$3\"" \
  -ex "detach" \
  -ex "quit"

# gdb -n -q -ex "attach $(pidof $1)" \ -ex "set $2 = \"$3\""


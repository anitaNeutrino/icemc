#!/usr/bin/env bash
target_src="$2"
linemarker="$3"
obj="$4"

lineno=$(grep -n "$linemarker" "$target_src" | gawk -F: '{print $1}' -)
echo $0:: Found GDB marker at line: "$lineno"

EXE="$1"
OUT="$5"
COND="$6"
rm "$OUT" # gdb would _append_ the log output otherwise.

gdb "$EXE" -q -batch \
  -ex "set print repeats 0" \
  -ex "set height 0" \
  -ex "set print elements 50000" \
  -ex "set logging file $OUT" \
  -ex "set logging on" \
  -ex "set logging redirect on" \
  -ex "b $target_src: $lineno $COND" \
  -ex "run" \
  -ex "echo BIZDOOMCAR\n" \
  -ex "output $obj"

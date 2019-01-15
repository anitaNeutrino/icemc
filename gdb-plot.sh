#!/usr/bin/env bash
gdb ./testEAS -q -batch \
  -ex "b ChanTrigger.cc:799" \
  -ex "run" \
  -ex "plot1d_opt tmp_volts[0] 'with lines'" \
  -ex "quit"

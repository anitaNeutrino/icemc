#!/usr/bin/env bash
gdb ./testEAS-1355 -q -batch \
  -ex "set print repeats 0" \
  -ex "set height 0" \
  -ex "set print elements 50000" \
  -ex "set logging file GDB/amp_eq_1/tmp_volts.dump" \
  -ex "set logging on" \
  -ex "set logging redirect on" \
  -ex "b ChanTrigger.cc:799" \
  -ex "run" \
  -ex "echo BIZDOOMCAR\n" \
  -ex "output tmp_volts[0]"

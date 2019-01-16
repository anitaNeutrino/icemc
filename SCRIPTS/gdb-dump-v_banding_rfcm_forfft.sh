#!/usr/bin/env bash
gdb ./testEAS-ae08 -q -batch \
  -ex "set print repeats 0" \
  -ex "set height 0" \
  -ex "set print elements 50000" \
  -ex "set logging file GDB/amp_eq_1/v_banding_rfcm_forfft.dump" \
  -ex "set logging on" \
  -ex "set logging redirect on" \
  -ex "b ChanTrigger.cc:938" \
  -ex "run" \
  -ex "echo BIZDOOMCAR\n" \
  -ex "output v_banding_rfcm_forfft[0][4]"

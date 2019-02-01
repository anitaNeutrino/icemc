#!/usr/bin/env bash

MD5_OUT_PREF=parent ~/PROGS/DepTrack/with-md5.sh --snapshot-only --inp $0
read SELF < val-parent-ID0-SNAPSHOT.out

EXE1=./testEAS

echo " ----------- ZhsTimeE -------- "
if true
then
EXE=SCRIPTS/gdb-dump-obj.sh
INP1=./.gdbinit

dump_name=ZhsTimeE.dump
dump_dir=GDB/DUMP
OUT="$dump_dir/$dump_name"

target_src=testEAS.cc
linemarker="GDB-DUMP-MARKER: ZhsTimeE"
obj=ZhsTimeE


MD5_OUT_PREF=zhs ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --devar HOME=$HOME \
 --devar ICEMC_SRC_DIR=$ICEMC_SRC_DIR \
 --devar LD_LIBRARY_PATH=$ANITA_UTIL_INSTALL_DIR/lib \
 --inp $EXE \
 --inp $INP1 \
 --inp $target_src \
 --out $OUT \
 "$EXE" "$EXE1" "$target_src" "$linemarker" "$obj" "$OUT"

read MD5 < val-zhs-ID0-MD5.out

EXE=SCRIPTS/column-from-parlist.sh
INP=$dump_dir/md5-${MD5}_$dump_name
OUT_BASE_NAME_COL1="ZhsTimeE-col.dat"
OUT=$dump_dir/$OUT_BASE_NAME_COL1

MD5_OUT_PREF=col ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --inp $EXE \
 --inp $INP \
 --out $OUT \
 $EXE $INP $OUT
fi
read MD5_COL1 < val-col-ID0-MD5.out

echo " ----------- volts_rx_forfft -------- "

EXE=SCRIPTS/gdb-dump-obj.sh
INP1=./.gdbinit

dump_name=vrx.dump
dump_dir=GDB/DUMP
OUT="$dump_dir/$dump_name"

target_src="ChanTrigger.cc"
linemarker="GDB-DUMP-MARKER: volts_rx_forfft"
obj="volts_rx_forfft[0][4][0]@512"
echo "$0: The obj: $obj"

COND="if ant == 30"
MD5_OUT_PREF=vrx ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --devar HOME=$HOME \
 --devar ICEMC_SRC_DIR=$ICEMC_SRC_DIR \
 --devar LD_LIBRARY_PATH=$ANITA_UTIL_INSTALL_DIR/lib \
 --inp $EXE \
 --inp $INP1 \
 --inp $target_src \
 --out $OUT \
 "$EXE" "$EXE1" "$target_src" "$linemarker" "$obj" "$OUT" "$COND"

read MD5 < val-vrx-ID0-MD5.out

JOB_NAME=vrx-col
OUT_BASE_NAME_COL2=vrx-col.dat
EXE=SCRIPTS/column-from-parlist.sh
INP=$dump_dir/md5-${MD5}_$dump_name
OUT=$dump_dir/$OUT_BASE_NAME_COL2

MD5_OUT_PREF=$JOB_NAME ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --inp $EXE \
 --inp $INP \
 --out $OUT \
 $EXE $INP $OUT

read MD5_COL2 < val-$JOB_NAME-ID0-MD5.out

echo " ----------- Comparison ------------- "

JOB_NAME=gpl
INP1=$dump_dir/md5-${MD5_COL1}_$OUT_BASE_NAME_COL1 # Using previous output name.
INP2=$dump_dir/md5-${MD5_COL2}_$OUT_BASE_NAME_COL2 # Using previous output name.
OUT_BASE_NAME=zhs-vs-vrx.eps
OUT=$dump_dir/$OUT_BASE_NAME
EXE=SCRIPTS/two-col-plot.sh

MD5_OUT_PREF=$JOB_NAME ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --inp $EXE \
 --inp $INP1 \
 --inp $INP2 \
 --out $OUT \
 $EXE $INP1 $INP2 $OUT "[500:1000]" "[]" "1e3" "+550" $dump_dir/_$JOB_NAME.gpl

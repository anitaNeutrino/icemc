#!/usr/bin/env bash

MD5_OUT_PREF=parent ~/PROGS/DepTrack/with-md5.sh --snapshot-only --inp $0
read SELF < val-parent-ID0-SNAPSHOT.out


EXE=SCRIPTS/gdb-dump-v_banding_rfcm_forfft.sh
INP=./.gdbinit
OUT=GDB/amp_eq_1/v_banding_rfcm_forfft.dump

MD5_OUT_PREF=rfcm ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --devar HOME=$HOME \
 --devar ICEMC_SRC_DIR=$ICEMC_SRC_DIR \
 --devar LD_LIBRARY_PATH=$ANITA_UTIL_INSTALL_DIR/lib \
 --inp $INP \
 --inp $EXE \
 --out $OUT \
 $EXE 

read MD5 < val-rfcm-ID0-MD5.out

EXE=SCRIPTS/column-from-parlist.sh
INP=GDB/amp_eq_1/md5-${MD5}_v_banding_rfcm_forfft.dump
OUT=GDB/amp_eq_1/v_banding_rfcm_forfft_col.dat

MD5_OUT_PREF=rfcm ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --inp $EXE \
 --inp $INP \
 --out $OUT \
 $EXE $INP $OUT

read MD5_RFCM_COL < val-rfcm-ID0-MD5.out
# exit 0


# read MD5 < val-rfcm-ID0-MD5.out
# INP=GDB/amp_eq_1/md5-${MD5}_v_banding_rfcm_forfft_col.dat
# OUT=GDB/amp_eq_1/tmp_volts-cmp.eps
# EXE=SCRIPTS/gdb-tmp_volts-plot.sh
# MD5_OUT_PREF=rfcm ~/PROGS/DepTrack/with-md5.sh \
#  --envar PARENT_SCRIPT=$SELF \
#  --envar PARENT_LINENO=$LINENO \
#  --devar PATH=/usr/bin \
#  --inp $EXE \
#  --inp $INP \
#  --out $OUT \
#  $EXE $INP $OUT


EXE=SCRIPTS/gdb-dump.sh

MD5_OUT_PREF=dump ~/PROGS/DepTrack/with-md5.sh \
 --devar PATH=/usr/bin \
 --devar HOME=$HOME \
 --devar ICEMC_SRC_DIR=$ICEMC_SRC_DIR \
 --devar LD_LIBRARY_PATH=$ANITA_UTIL_INSTALL_DIR/lib \
 --inp $EXE \
 --out GDB/amp_eq_1/tmp_volts.dump \
 $EXE 

read MD5 < val-dump-ID0-MD5.out
# echo "MD5 is: $MD5"
# ls -lh GDB/amp_eq_1/md5-${MD5}_tmp_volts.dump

EXE=SCRIPTS/column-from-parlist.sh
INP=GDB/amp_eq_1/md5-${MD5}_tmp_volts.dump
OUT=GDB/amp_eq_1/tmp_volts_col.dat

MD5_OUT_PREF=dump ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --inp $EXE \
 --inp $INP \
 --out $OUT \
 $EXE $INP $OUT


read MD5 < val-dump-ID0-MD5.out
INP1=GDB/amp_eq_1/md5-${MD5}_tmp_volts_col.dat
INP2=GDB/amp_eq_1/md5-${MD5_RFCM_COL}_v_banding_rfcm_forfft_col.dat
OUT=GDB/amp_eq_1/tmp_volts-cmp.eps
EXE=SCRIPTS/gdb-tmp_volts-plot.sh
MD5_OUT_PREF=dump ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --inp $EXE \
 --inp $INP1 \
 --inp $INP2 \
 --out $OUT \
 $EXE $INP1 $INP2 $OUT


read MD5 < val-dump-ID0-MD5.out
INP=GDB/amp_eq_1/md5-${MD5}_tmp_volts-cmp.eps
SCRIPTS/view-gpl-eps.sh $INP

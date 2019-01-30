#!/usr/bin/env bash

MD5_OUT_PREF=parent ~/PROGS/DepTrack/with-md5.sh --snapshot-only --inp $0
read SELF < val-parent-ID0-SNAPSHOT.out

EXE=SCRIPTS/testEAS-save-std.sh
# EXE1=./testEAS-ec15
# EXE1=./testEAS-6e59
# EXE1=./testEAS-c04f
# EXE1=./testEAS-85b5
# EXE1=./testEAS-c929
# EXE1=./testEAS-e2a1
# EXE1=./testEAS-1438
EXE1=./testEAS-5fe3
INP1=inputs.anita3.conf
INP2=/nfs/data_disks/herc0a/users/bugaev/ANITA/SIMS/Event_4212/timefresnel-root.dat
OUT=PLOTS/WForms/MD5/testEAS-std.txt

# MD5_OUT_PREF=std ~/PROGS/DepTrack/with-md5.sh --dry-run-out dry
MD5_OUT_PREF=std ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --devar HOME=$HOME \
 --devar ICEMC_SRC_DIR=$ICEMC_SRC_DIR \
 --devar LD_LIBRARY_PATH=$ANITA_UTIL_INSTALL_DIR/lib \
 --inp $EXE \
 --inp $INP1 \
 --inp $INP2 \
 --out $OUT \
 $EXE $EXE1 $INP1 $INP2 $OUT

# exit 0

read MD5 < val-std-ID0-MD5.out

EXE=grep
INP=PLOTS/WForms/MD5/md5-${MD5}_testEAS-std.txt
OUT=PLOTS/WForms/MD5/testEAS-std-only-hv.txt

MD5_OUT_PREF=std ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --inp $INP \
 --out $OUT \
 bash -c "$EXE Antenna $INP > $OUT"
 
read MD5 < val-std-ID0-MD5.out

EXE=SCRIPTS/transpose.sh
INP=PLOTS/WForms/MD5/md5-${MD5}_testEAS-std-only-hv.txt
OUT=PLOTS/WForms/MD5/testEAS-hv-cols.txt
MD5_OUT_PREF=std ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --inp $EXE \
 --inp $INP \
 --out $OUT \
bash -c "$EXE $INP > $OUT"

read MD5 < val-std-ID0-MD5.out

EXE=SCRIPTS/group-prefixed-rows.sh
INP=PLOTS/WForms/MD5/md5-${MD5}_testEAS-hv-cols.txt
OUT=PLOTS/WForms/MD5/h.txt
MD5_OUT_PREF=hor ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --inp $EXE \
 --inp $INP \
 --out $OUT \
$EXE $INP 'H:' $OUT

EXE=SCRIPTS/group-prefixed-rows.sh
INP=PLOTS/WForms/MD5/md5-${MD5}_testEAS-hv-cols.txt
OUT=PLOTS/WForms/MD5/v.txt
MD5_OUT_PREF=ver ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --inp $EXE \
 --inp $INP \
 --out $OUT \
$EXE $INP 'V:' $OUT

read MD5HOR < val-hor-ID0-MD5.out
read MD5VER < val-ver-ID0-MD5.out

EXE=/nfs/data_disks/herc0a/users/bugaev/ANITA/anitaBuildTool/components/icemc/PLOTS/WForms/H-V-WF.sh
INP1=PLOTS/WForms/MD5/md5-${MD5HOR}_h.txt
INP2=PLOTS/WForms/MD5/md5-${MD5VER}_v.txt
TMPSCRIPT=PLOTS/WForms/MD5/_all-antennas.gpl
OUT=PLOTS/WForms/MD5/H-V-WF.eps
MD5_OUT_PREF=gpl ~/PROGS/DepTrack/with-md5.sh \
 --envar PARENT_SCRIPT=$SELF \
 --envar PARENT_LINENO=$LINENO \
 --devar PATH=/usr/bin \
 --inp $EXE \
 --inp $INP1 \
 --inp $INP2 \
 --out $OUT \
$EXE $INP1 $INP2 $TMPSCRIPT $OUT

#Comment.
#More nonsense.

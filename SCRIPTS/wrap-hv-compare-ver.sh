#!/usr/bin/env bash

MD5_OUT_PREF=parent ~/PROGS/DepTrack/with-md5.sh --snapshot-only --inp $0
read SELF < val-parent-ID0-SNAPSHOT.out


MD5OLD=50d47adc9b0a25e71270de28f050c5db # H-pol in ANITA coll call.
# MD5NEW=2f49e8b9150977787598e6ac77c9b9f7 # H-pol after propagation of phases.
MD5NEW=78e7079616bf8a0febf5ba51302b0632 # H-pol after propagation of phases, fully recompiled.

EXE=/nfs/data_disks/herc0a/users/bugaev/ANITA/anitaBuildTool/components/icemc/SCRIPTS/H-H-WF.sh
INP1=PLOTS/WForms/MD5/md5-${MD5OLD}_h.txt
INP2=PLOTS/WForms/MD5/md5-${MD5NEW}_h.txt
TMPSCRIPT=PLOTS/WForms/MD5/_all-antennas.gpl
OUT=PLOTS/WForms/MD5/H-H-WF.eps
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

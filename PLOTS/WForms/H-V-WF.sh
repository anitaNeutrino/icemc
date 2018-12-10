INP1=$1
INP2=$2
TMP_SCRIPT=$3
OUT=$4

# TMP_SCRIPT: PLOTS/WForms/MD5/_all-antennas.gpl
cat << EOF > $TMP_SCRIPT 
set out "$OUT"
set term post color solid 5 

set multiplot layout 6,8
set macro
do for [i = 1:48] {
   # set nokey
   if (i == 31) {
     nplot = sprintf("Selected: ------> %d <------", i)
   }
   else {
     nplot = sprintf("%d", i)
   }
   set title nplot
 plot [][-3e+1 : +3e+1] \
 "$INP1" using 0:1  with lines title "H",\
 "$INP2" using 0:1  with lines title "V"
 # "$INP1" using 0:(column(i) / 1e6)  with lines title "H",\
 # "$INP2" using 0:(column(i) / 1e6)  with lines title "V"
  
}
unset multiplot
EOF

gnuplot "$TMP_SCRIPT"

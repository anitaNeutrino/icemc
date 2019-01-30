INP1=$1
INP2=$2
TMP_SCRIPT=$3
OUT=$4

# TMP_SCRIPT: PLOTS/WForms/MD5/_all-antennas.gpl
cat << EOF > $TMP_SCRIPT 
set out "$OUT"
set term post color solid 5 

# set multiplot layout 6,8
set multiplot layout 1,1
set macro
# do for [i = 1:48] {

do for [i = 31:31] {
   # set nokey
   if (i == 31) {
     nplot = sprintf("Selected: ------> %d <------", i)
   }
   else {
     nplot = sprintf("%d", i)
   }
   set title nplot

# plot [][-100 : +100]\

plot [][]\
	 "$INP1" using 0:(column(i))  with lines title "phase: f",\
	 "$INP2" using 0:(column(i))  with lines title "phase: t"
  
}
unset multiplot
EOF

gnuplot "$TMP_SCRIPT"

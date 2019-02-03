#!/usr/bin/env bash

NINP=0
INP=none
OUT=none
XRANGE=none
YRANGE=none
TMP_GNUPLOT_SCRIPT=none

function display_help()
{
    echo "Usage: $0 --option1 --option2" 1>&2
}


if [[ $# == 0 ]]
then
    echo "Expected input parameters are missing." 1>&2
    display_help
    exit 100 # Mandatory input parameters are missing.
fi

function reset()
{
    XSCALE=1.0
    YSCALE=1.0
    XSHIFT=+0.0
    YSHIFT=+0.0
    LINETYPE=none
    LINEWIDTH=none
}

PLOTACC="plot"
LINEACC=""

function flush_line()
{
    local LINE_SEP="$1"
    local LT=""
    test "$LINETYPE" != "none" && LT=" lt $LINETYPE"
    local LW=""
    test "$LINEWIDTH" != "none" && LW=" lw $LINEWIDTH"
    local APPENDIX="\"$INP\" (column(0) * $XSCALE $XSHIFT):(column(1) * $YSCALE $YSHIFT) with line$LT$LW""$LINE_SEP"
    LINEACC="${LINEACC}$APPENDIX" 
}

reset

while [[ $# > 0 ]]
do
    case "$1" in
      --xrange)
         XRANGE="$2"
         shift 2
         ;;
      --yrange)
         YRANGE="$2"
         shift 2
         ;;
      --lt)
         LINETYPE="$2"
         shift 2
         ;;
      --lw)
         LINEWIDTH="$2"
         shift 2
         ;;
      --xshift)
         XSHIFT="$2"
         shift 2
         ;;
      --yshift)
         YSHIFT="$2"
         shift 2
         ;;
      --xscale)
         XSCALE="$2"
         shift 2
         ;;
      --yscale)
         YSCALE="$2"
         shift 2
         ;;
      --linkify)
        LINKIFY=yes
        shift 1
        ;;
      --inp) # Input file.
         NINP=$((NINP + 1))
         if [[ "$NINP" > 1 ]] # We are done with the previous input file, flush its plotting line.
         then
             LINE_SEP=$',\\\n'
             flush_line "$LINE_SEP"
             reset
         fi
         INP="$2"
         shift 2
         ;;
      --out) # Output file.
         OUT="$2"
         shift 2
         ;;
      -h | --help)
	  display_help  # Call your function
	  # no shifting needed here, we're done.
	  exit 0
	  ;;
      -v | --verbose)
	  verbose="verbose"
	  shift
	  ;;
      --) # End of all options
	  shift
	  break;;
      -*)
	  echo "Error: Unknown option: $1" >&2
          display_help
	  exit 107 
	  ;;
      *)  # No more options
	  break
	  ;;
    esac
done

LINE_SEP="" # Empty
flush_line "$LINE_SEP"

YR=""
test "$YRANGE" != "none" && YR=" $YRANGE"
XR=""
test "$XRANGE" != "none" && XR=" $XRANGE"
test "$YRANGE" != "none" && test "$XRANGE" = "none" && XR=" []"
PLOTACC="$PLOTACC$XR$YR"

ACC="$PLOTACC "$'\\\n'"$LINEACC"

echo "$ACC"

exit 0

echo "Temporary gnuplot script is: $TMP_GNUPLOT_SCRIPT"

cat << EOF  > "$TMP_GNUPLOT_SCRIPT"
set out "$OUT"
set term post color solid 18
plot $XRANGE $YRANGE \
       	"$FILE1" using 0:(column(1) * 1.0) with lines  lt 1 lw 3,\
       	"$FILE2" using (column(0) $SHIFT):(column(1) * $SCALE) with lines lt 2
EOF

gnuplot $TMP_GNUPLOT_SCRIPT

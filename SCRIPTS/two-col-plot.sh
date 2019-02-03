#!/usr/bin/env bash

FILE1=$1
FILE2=$2
OUT=$3
XRANGE=$4
YRANGE=$5
SCALE=$6
SHIFT=$7
TMP_GNUPLOT_SCRIPT=$8

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

while [[ $# > 0 ]]
do
    case "$1" in
      --linkify)
        LINKIFY=yes
        shift 1
        ;;
      --inp) # Input file.
         # Either a dependency file or a file with a list of files to be pointed by symbolic links in the directory specified by the --out parameter.

         INPS+=("$2")
         INPCNT=$((INPCNT + 1))
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

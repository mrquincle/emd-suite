#!/bin/sh

usage="$0 <input directory> <output directory>"

idir=${1:? "Usage: $usage"}
odir=${2:? "Usage: $usage"}

mkdir -p $odir

# iterate over each directory
for D in $idir/*/; do
  #echo $D
  parent=$(basename $D)
  latest=$(echo $idir/$parent/LATEST)
  rfile=$(echo $parent.txt)
  ./organize_for_plot.sh $latest $rfile
  ./display_cube.R $latest/$rfile $odir
done

mkdir -p $odir/summary
pdftk $odir/*.pdf cat output $odir/summary/results.pdf

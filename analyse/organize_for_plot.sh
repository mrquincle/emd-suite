#!/bin/bash

usage="$0 <directory>"

dir=${1:? "Usage: $usage"}

ofile=${2:-results_for_plot.txt}

shopt -s extglob

cd $dir

rm -f $ofile

for f in results+([0-9]).txt; do
  #echo $f
  [[ "$f" =~ ([0-9]+) ]] && value="${BASH_REMATCH[1]}"
  #echo $value
 
  cat $f | tr ' ' '\t' | awk -v value="$value" -F "\t" '{print $0,value}' >> $ofile
done

echo "Stored $(wc -l $ofile | cut -f1 -d' ') items in $ofile"

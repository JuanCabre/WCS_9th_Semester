#!/bin/bash
#Convert images to eps.
#Usage: ./fig_conv.sh (figure to convert) (name of the new figure no extension)

output=$2
output+=".pdf"

echo $output
convert $1 $output
pdftops -eps $output

rm $output

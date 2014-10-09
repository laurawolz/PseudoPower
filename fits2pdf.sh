#!/bin/sh
echo -e "\033[34m Type './fits2pdf help' if you need help \033[0m"
if [ "$1" == "help" ] ;
then echo -e "\033[35m First input parameter must be path+filename without . and ending; second and third input is min and max value of color palette (optional); \033[0m "
else
INFILE=$1
MIN=$2
MAX=$3

echo $INFILE.fits to $INFILE.pdf
echo start
 map2tga $INFILE.fits $INFILE.tga -bar  # -max $MAX -min $MIN -pal 4 # -log  # -add 40 -log
 #map2tga $INFILE.fits $INFILE.tga  -bar

convert $INFILE.tga $INFILE.pdf
rm $INFILE.tga
fi

#!/bin/bash

OUTFILE=$(echo $1 | cut -f 1 -d '.')
GENOME=$2

cat ${GENOME}_header.txt > $OUTFILE.pairs
echo "#command: bash bedpe2cool.sh" "$@" >> $OUTFILE.pairs
# convert to pairs file 
python bedpe2pairs.py $1 $OUTFILE.pairs

# sort and zip output file
sort -k2,2 -k4,4 -k3,3n -k5,5n $OUTFILE.pairs | bgzip -c > $OUTFILE.pairs.gz
# index zipped file
pairix -f $OUTFILE.pairs.gz

# generate cool file
cooler cload pairix $GENOME.chrom.sizes.mainonly:1000000000 $OUTFILE.pairs.gz $OUTFILE.gini.cool
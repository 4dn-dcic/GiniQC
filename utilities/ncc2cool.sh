#!/bin/bash
:<<'end_long_comment'
---- COPYRIGHT ----------------------------------------------------------------
Copyright (C) 2020
Connor Horton (Harvard University)

---- LICENSE ------------------------------------------------------------------
This file is part of GiniQC.

GiniQC is free software: you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

GiniQC is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details. 

You should have received a copy of the GNU Lesser General Public License along
with this software.  If not, see <http://www.gnu.org/licenses/>.
end_long_comment


OUTFILE=$(echo $1 | cut -f 1 -d '.')
GENOME=$2

cat ${GENOME}_header.txt > $OUTFILE.pairs
echo "#command: bash ncc2cool.sh" "$@" >> $OUTFILE.pairs
# convert to pairs file 
python ncc2pairs.py $1 $OUTFILE.pairs

# sort and zip output file
sort -k2,2 -k4,4 -k3,3n -k5,5n $OUTFILE.pairs | bgzip -c > $OUTFILE.pairs.gz
# index zipped file
pairix -f $OUTFILE.pairs.gz

# generate cool file
cooler cload pairix $GENOME.chrom.sizes.mainonly:1000000000 $OUTFILE.pairs.gz $OUTFILE.gini.cool
#!/bin/bash

:<<'end_long_comment'
---- COPYRIGHT ----------------------------------------------------------------
Copyright (C) 2017-2018
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

usage() {
	cat <<EOM
Usage: bash GiniQC.sh [-h] -f FILE(s) -o OUTFILE -b BEDFILE [-c CISTHRESHOLD] [-g GINITHRESHOLD] [-r READSTHRESHOLD] [-a MAXABERRATION]
	-h|--help				prints this message
	-f FILE(s)			path to cooler matrix file (must end in .cool) or path to a list of files (any other extension)
	-o OUTFILE			desired name for output files, including pairs file
	-b BEDFILE			paired-end fastq files corresponding to a single cell
	-c CISTHRESHOLD		user-defined value (default: 80% cis)
	-g GINITHRESHOLD	minimum GiniQC value (if unspecified, our tool will suggest a threshold)
	-r READSTHRESHOLD 	minimum number of reads per cell (default: 10,000 reads)
	-a MAXABERRATION 	maximum fold-change in coverage  (default: 2-fold)
EOM
	exit 1;
}
if [ $# == 0 ] ; then
    usage
    exit 1;
fi

CISTHRESHOLD=80
READSTHRESHOLD=10000
MAXABERRATION=2

while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -h|--help)
	    usage
	    shift # past argument
	    shift # past value
	    ;;
    -o)
	    OUTFILE="$2"
	    shift # past argument
	    shift # past value
	    ;;
    -f)
	    INFILE="$2"
	    shift # past argument
	    shift # past value
	    ;;
    -b)
		CHROMS="$2"
        shift
        shift
        ;;
    -c)
		CISTHRESHOLD="$2"
        shift
        shift
        ;;
    -g)
		GINITHRESHOLD="$2"
        shift
        shift
        ;;
    -r)
		READSTHRESHOLD="$2"
        shift
        shift
        ;;
    -a)
		MAXABERRATION="$2"
        shift
        shift
        ;;
    *)    # unknown option
    	POSITIONAL+=("$1") # save it in an array for later
    	shift # past argument
    	;;
esac
done

# if user provides a single .cool file, output statistics to single text file
if [[ $INFILE =~ ".cool"$ ]] ; then
	echo "running GiniQC on $FILE"
	python gini.py $INFILE $OUTFILE
# if user provides a list of .cool files and no predetermined Gini threshold
elif [ -z "$GINITHRESHOLD" ] ;
	echo "determining a threshold for GiniQC"
	python threshold.py $INFILE $CHROMS $OUTFILE > temp.out
	GINITHRESHOLD=$(cat temp.out)
	rm temp.out
	echo "running GiniQC for each file in $FILE"
	python gini_wrapper.py $INFILE $OUTFILE $READSTHRESHOLD $CISTHRESHOLD $MAXABERRATION $GINITHRESHOLD
# if user provides a list of .cool files and a Gini threshold
else
	echo "running GiniQC for each file in $FILE"
	python gini_wrapper.py $INFILE $OUTFILE $READSTHRESHOLD $CISTHRESHOLD $MAXABERRATION $GINITHRESHOLD
fi
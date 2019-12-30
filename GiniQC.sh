#!/bin/bash

:<<'end_long_comment'
---- COPYRIGHT ----------------------------------------------------------------
Copyright (C) 2017-2020
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
Usage: bash GiniQC.sh [-h] -f FILE(s) -o OUTFILE -b BEDFILE [-c CISTHRESHOLD] [-g GINITHRESHOLD] [-r READSTHRESHOLD] [-a MAXABERRATION] [-x]
	-h|--help				prints this message
	-f FILE(s)			path to cooler matrix file (must end in .cool) or path to a list of files (any other extension)
	-o OUTFILE			desired name for output files
	-b BEDFILE			paired-end fastq files corresponding to a single cell
	-c CISTHRESHOLD		minimum percent cis value per cell (default: 80)
	-g GINITHRESHOLD	minimum GiniQC value (if not specified, our tool will suggest a threshold)
	-r READSTHRESHOLD 	minimum number of reads per cell (default: 10,000 reads)
	-a MAXABERRATION 	maximum fold-change in chromosomal sequencing coverage (default: 2-fold)
	-x					when used, GiniQC threshold is determined only on cells passing cis threshold (see -c above)
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
CULL_BY_CIS=false

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
    -x)
		CULL_BY_CIS=true
		shift
		shift
		;;
    *)    # unknown option
    	POSITIONAL+=("$1") # save it in an array for later
    	shift # past argument
    	;;
esac
done

if [[ $INFILE =~ ".cool"$ ]] ; then # if user provides a single .cool file, output statistics to single text file
	echo ">> Running GiniQC on ${INFILE}"
	python gini.py $INFILE $OUTFILE
elif [ -z "$GINITHRESHOLD" ] ; then # if user provides a list of .cool files and no predetermined Gini threshold
	echo ">> Determining a threshold for GiniQC"
	if [[ "$CULL_BY_CIS" == "true" ]]; then
		GINITHRESHOLD=$(python threshold.py $INFILE $CHROMS $CISTHRESHOLD)
	else
		GINITHRESHOLD=$(python threshold.py $INFILE $CHROMS)
	fi
	echo ">> Running GiniQC for each file in ${INFILE}"
	python gini_wrapper.py $INFILE $OUTFILE $READSTHRESHOLD $CISTHRESHOLD $MAXABERRATION $GINITHRESHOLD
else # if user provides a list of .cool files and a Gini threshold
	echo ">> Running GiniQC for each file in ${INFILE}"
	python gini_wrapper.py $INFILE $OUTFILE $READSTHRESHOLD $CISTHRESHOLD $MAXABERRATION $GINITHRESHOLD
fi

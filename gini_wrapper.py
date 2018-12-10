"""
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
"""

import os, sys
import cooler
from gini import normalize_matrix, gini, adjust, get_max_aberration

OUTFILE_HEADER = "file_name\treads\t% cis\tGiniQC(raw)\tGiniQC(adj)\tMax coverage aberration\tPassed?\n"
LINE_STRUCTURE = "%s\t%d\t%.2f\t%.3f\t%.3f\t%f\t%s\n"

def main():
	# parse arguments
	infile = open(sys.argv[1], 'r')
	outfile = open(sys.argv[2], 'w')
	reads_threshold = int(sys.argv[3])
	cis_threshold = float(sys.argv[4])
	max_aberration = float(sys.argv[5])
	gini_threshold = float(sys.argv[6].strip())

	# print header lines
	outfile.write(OUTFILE_HEADER)
	outfile.write(LINE_STRUCTURE % ("THRESHOLD", reads_threshold, cis_threshold, 0.0, gini_threshold, max_aberration, "N/A"))

	# calculate GiniQC and other metrics for each cell. Then write each to outfile
	for line in infile:
		file = line.strip()
		cell_name = line.split(".cool")[0]
		matrix = cooler.Cooler(file)
		normalized, reads, cis_reads, trans_reads = normalize_matrix(matrix)
		percent_cis = 100.0*cis_reads/reads
		raw_gini = gini(normalized)
		adj_gini = adjust(raw_gini, reads)
		chrom_aberration = get_max_aberration(matrix)

		passed = (reads > reads_threshold) and (percent_cis > cis_threshold) and (adj_gini > gini_threshold) and (chrom_aberration < max_aberration)
		outfile.write(LINE_STRUCTURE % (cell_name, reads, percent_cis, raw_gini, adj_gini, chrom_aberration, passed))


if __name__ == '__main__':
	main()
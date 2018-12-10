import os, sys
import cooler
from gini import normalize_matrix, gini, adjust, get_max_aberration

OUTFILE_HEADER = "file_name\tnum_reads\t%_cis\tGiniQC(raw)\tGiniQC(adjusted)\tMax_coverage_aberration\tPassed?\n"
LINE_STRUCTURE = "%s\t%d\t%.2f\t%.3f\t%.3f\t%f\t%s"

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
	outfile.write(LINE_STRUCTURE % ("TRESHOLD", reads_threshold, cis_threshold, 0.0, gini_threshold, max_aberration, "N/A"))

	# calculate GiniQC and other metrics for each cell. Then write each to outfile
	for line in infile:
		file = line.strip()
		cell_name = line.split(".cool")[0]
		matrix = cooler.Cooler(file)
		normalized, reads, cis_reads, trans_reads = normalize_matrix(matrix)
		cis = cis_reads/reads
		raw_gini = gini(normalized)
		adj_gini = adjust(raw_gini, reads)
		max_aberration = get_max_aberration(matrix)

		passed = (reads > reads_threshold) and (cis > cis_threshold) and (adj_gini > gini_threshold) and (chrom_aberration < max_aberration)
		outfile.write(LINE_STRUCTURE % (cell_name, reads, cis, raw_gini, adj_gini, chrom_aberration, passed))


if __name__ == '__main__':
	main()
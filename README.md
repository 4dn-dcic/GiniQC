# GiniQC

#### Dependencies

- Python 3
- Cooler (`pip install cooler`)
- numpy (`pip install numpy`)
- pandas (`pip install pandas`)

#### Usage

```Usage: bash GiniQC.sh [-h] -f FILE(s) -o OUTFILE -b BEDFILE [-c CISTHRESHOLD] [-g GINITHRESHOLD] [-r READSTHRESHOLD] [-a MAXABERRATION]
	-h|--help				prints this message
	-f FILE(s)			path to cooler matrix file (must end in .cool) or path to a list of files (any other extension)
	-o OUTFILE			desired name for output files, including pairs file
	-b BEDFILE			paired-end fastq files corresponding to a single cell
	-c CISTHRESHOLD		user-defined value (default: 80% cis)
	-g GINITHRESHOLD	minimum GiniQC value (if unspecified, our tool will suggest a threshold)
	-r READSTHRESHOLD 	minimum number of reads per cell (default: 10,000 reads)
	-a MAXABERRATION 	maximum fold-change in coverage  (default: 2-fold)```
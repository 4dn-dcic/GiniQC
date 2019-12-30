#pairs columns: readID chr1 pos1 chr2 pos2 strand1 strand2
#bedpe columns: chr1 start1 end1 chr2 start2 end2 name score strand1 strand2
import sys

def main():
	infile = open(sys.argv[1], 'r')
	outfile = open(sys.argv[2], 'a')

	for i, line in enumerate(infile):
		chr_a, start_a, end_a, chr_b, start_b, end_b, name, score, strand_a, strand_b = line.strip().split('\t')[:10]
		out = [name]
		if strand_a == "+":
			out.extend([chr_a, start_a])
		else:
			out.extend([chr_a, end_a])
		if strand_b == "+":
			out.extend([chr_b, start_b])
		else:
			out.extend([chr_b, end_b])
		out.extend([strand_a, strand_b])
		outfile.write("\t".join(out)+"\n")

	infile.close()
	outfile.close()


if __name__ == '__main__':
	main()
import numpy as np
import pandas as pd
import os, sys, cooler
from itertools import combinations
from collections import defaultdict

from gini_zeros import normalize_matrix, gini, adjust

def make_df(bedfile):
	chroms, starts, ends = [], [], []
	for line in bedfile:
		chrom, start, end = line.strip('\n').split('\t')
		chroms.append(chrom)
		starts.append(start)
		ends.append(end)
	return pd.DataFrame({'chrom':pd.Series(chroms, dtype='str'),
		'start':pd.Series(starts, dtype='int'),
		'end':pd.Series(ends, dtype='int')})


def main():
	filelist = open(sys.argv[1], 'r')
	bedfile = open(sys.argv[2])
	bins_df = make_df(bedfile)
	files = []
	for line in filelist:
		files.append(line.strip('\n'))
	numfiles = len(files)

	np.random.seed(0)

	if numfiles < 40:
		combos = list(combinations(files, 2))
		if len(combos) > 20:
		 	combos = [combos[x] for x in np.random.choice(len(combos), 20, replace=False)]

	else:
		selection = np.random.choice(len(files), 400, replace=True)
		combos = np.reshape([files[x] for x in selection], (200,2))

	reads = {}
	rawgini = {}
	adjustedgini = {}

	for pair in combos:
		if pair[0] == pair[1]:
			continue
		pair = tuple(pair)
		try:
			cool1 = cooler.Cooler(pair[0]+".chrom.cool")
			cool2 = cooler.Cooler(pair[1]+".chrom.cool")
			matrix1 = np.array(cool1.matrix(as_pixels=True, balance=False)[:])
			matrix2 = np.array(cool2.matrix(as_pixels=True, balance=False)[:])
		except:
			continue
		numreads1 = sum(matrix1[:,-1])
		numreads2 = sum(matrix2[:,-1])

		if numreads1 == 0 or numreads2 == 0:
			continue

		totalreads = numreads1 + numreads2

		numtoselect = int(abs(np.random.normal(totalreads/2, totalreads/20)))

		rands = np.random.choice(np.arange(1,totalreads), numtoselect, replace=False)
		rands.sort()
		mixed_matrix = defaultdict(int)
		i = 0
		j = matrix1[i,-1]
		for rand in rands:
			if rand <= numreads1:
				while rand > j:
					i += 1
					j += matrix1[i,-1]
				else:
					bin1, bin2, count = matrix1[i]
					mixed_matrix[(bin1, bin2)] += 1
			else:
				while rand > j:
					i += 1
					j += matrix2[i-len(matrix1)-1,-1]
				else:
					bin1, bin2, count = matrix2[i-len(matrix1)]
					mixed_matrix[(bin1, bin2)] += 1

		mixedmatrix_columndict = {'bin1_id': [], 'bin2_id': [], 'count': []}
		for key, value in mixed_matrix.items():
			mixedmatrix_columndict['bin1_id'].append(key[0])
			mixedmatrix_columndict['bin2_id'].append(key[1])
			mixedmatrix_columndict['count'].append(value)
		
		pixel_df = pd.DataFrame(mixedmatrix_columndict).sort_values(by=['bin1_id', 'bin2_id']).reset_index(drop=True)
		cooler.io.create("temp.cool", bins=bins_df, pixels=pixel_df, dtype={'bin1_id':int, 'bin2_id':int, 'count':int})
		newcool = cooler.Cooler("temp.cool")

		normalized, reads[pair], cis, trans = normalize_matrix(newcool)
		rawgini[pair] = gini(normalized)
		adjustedgini[pair] = adjust(rawgini[pair], reads[pair])
		os.unlink("temp.cool")

	outfile = open("out.txt", 'w')
	outfile.write(",".join([str(v) for v in adjustedgini.values()]))
	outfile.close()


if __name__ == '__main__':
	main()
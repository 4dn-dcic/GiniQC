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

SLOPE = -0.0579
INTERCEPT = 1.233

import numpy as np
import sys, cooler, math
from collections import defaultdict

def get_max_abberation(cool):
    matrix_data = np.array(cool.matrix(as_pixels=True, balance=False, join=True)[:])
    matrix_bins = np.array(cool.bins()[:])

    chromsizes = {}
    for bin_ in matrix_bins:
        chrom, start, end, weight = bin_
        chromsizes[chrom] = end+1

    chrom_counts = defaultdict(int)
    for line in matrix_data:
        chr1, start1, end1, chr2, start2, end2, count = line
        if chr1 == chr2:
            chrom_counts[chr1] += count
        else:
            chrom_counts[chr1] += count
            chrom_counts[chr2] += count

    for chrom in chrom_counts:
        chrom_counts[chrom] /= chromsizes[chrom]
    coverage = list(chrom_counts.values())
    median_coverage = np.median(coverage)

    max_abberation = 1
    for chrom in coverage:
        if chrom > median_coverage:
            max_abberation = max(max_abberation, chrom/median_coverage)
        else:
            max_abberation = max(max_abberation, median_coverage/chrom)

    return max_abberation

def gini(values):
    values.sort()
    numerator, denominator = 0, 0
    count = 1
    n = len(values)
    for i in values:
        numerator += (2*count-n-1)*i
        denominator += i
        count += 1
    denominator *= n
    return numerator/denominator

def adjust(gini_index, numreads):
    if gini_index > 1e-4:
        return math.atan(gini_index - SLOPE*np.log(numreads) - INTERCEPT)/math.pi + 0.5
    else:
        return np.nan

def normalize_matrix(cool):
    # normalizes matrix by iterative correction
    bias, stats = cooler.ice.iterative_correction(cool, store=True, max_iters=200)
    # generates an array containing number of contacts, normalized weights for upper triangle of matrix
    matrix_data = np.array(cool.matrix(as_pixels=True, join=True, balance=True)[:])
    # gives total number of reads 
    total_sum = sum(matrix_data[:,-2])

    bins = {}
    count = 0

    nbins = cool.bins().shape[0]
    normalized = np.zeros((nbins, nbins), dtype=np.float64)

    last_chrom = 0
    first_chrom_bin = 0

    for bin_ in np.array(cool.bins()[:]):
        chrom, start, end, weight = bin_
        if chrom != last_chrom:
            # mark all of the cis bins as -1 in order to separate them out later
            for i in range(first_chrom_bin, count):
                for j in range(first_chrom_bin, count):
                    normalized[i,j] = -1
            # reset values for new chromosome
            first_chrom_bin = count
            last_chrom = chrom
        bins[(chrom, start, end)] = count
        count += 1
    # make sure cis bins are set at -1 for last chromosome
    for i in range(first_chrom_bin, count):
        for j in range(first_chrom_bin, count):
            normalized[i,j] = -1

    # select trans reads only
    cistotal = 0
    transtotal = 0
    for matrix_bin in matrix_data:
        chr1, start1, end1, chr2, start2, end2, count, weight = matrix_bin
        bin1 = bins[(chr1, start1, end1)]
        bin2 = bins[(chr2, start2, end2)]
        if chr1 != chr2:
            normalized[bin1, bin2] = weight
            transtotal += count
        else:
            cistotal += count

    #remove cis contacts, bottom diagonal before gini index calculation
    gini_array = []
    for i in range(nbins):
        for j in range(i+1, nbins):
            if normalized[i,j] >= 0:
                gini_array.append(normalized[i,j])
    return gini_array, total_sum, cistotal, transtotal

def main():
    matrix = cooler.Cooler(sys.argv[1]) # .cool filehandle
    # find maximum chromosomal read count aberration
    max_abberation = get_max_abberation(matrix)
    # first normalize matrix
    normalized, total_reads, cis_reads, trans_reads = normalize_matrix(matrix)
    # then calculate gini index
    if total_reads < 1:
        print("Input matrix %s is empty" % sys.argv[1])
    else:
        gini_idx = gini(normalized)
        adjusted_gini = adjust(gini_idx, total_reads)
        print("Gini coefficient for %s: %.3f" % (sys.argv[1], gini_idx))
        print("Adjusted Gini for %s: %.3f" % (sys.argv[1], adjusted_gini))
        if len(sys.argv) > 2:
            summary = open(sys.argv[2], 'r') # if used in conjunction with PairsQC or pipeline
            file_contents = summary.readlines()
            summary.close()
            # print Gini to summary text file created by PairsQC
            out = open(sys.argv[2], 'w')
            for line in file_contents:
                if "Gini" not in line:
                    out.write(line)
            out.write("%s\t%.3f\n" % ("Gini coefficient", gini_idx))
            out.write("%s\t%.3f\n" % ("Adjusted Gini coefficient", adjusted_gini))
            out.close()

if __name__ == '__main__':
    main()
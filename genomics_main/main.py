import pysam
import os
import numpy as np
import sys

outfile = open('outfile.txt', 'w')

vcf = pysam.VariantFile('experiment_0.vcf')
samples = np.array(vcf.header.samples)

# Question 1:
print('\nQuestion 1', file=outfile)
print('Number of people in the experiment: {}'.format(len(samples)), file=outfile)
infected_num = 0
for sample in samples:
    if sample.startswith('case'):
        infected_num += 1
print('Number of infected people: {}'.format(infected_num), file=outfile)
print('Number of healthy people: {}'.format(len(samples)-infected_num), file=outfile)

# Question 2:
print('\nQuestion 2',file=outfile)
variants = list(vcf.fetch())
print ('Number of mutations in total:', len(variants), file=outfile)
chrVariants = {}
for v in variants:
    if chrVariants.get(v.chrom) is None:
        chrVariants[v.chrom] = 0
    chrVariants[v.chrom]+=1
for chrm in chrVariants.keys():
    print('chromosome {} | mutations num {}'.format(chrm, chrVariants[chrm]),file=outfile)

outfile.close()

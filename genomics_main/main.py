import os
import pysam
import numpy as np
from matplotlib import pyplot as plt

currentDirectory = os.path.dirname(os.path.realpath(__file__))
print(currentDirectory)
outfile = open(currentDirectory + '/outfile.txt', 'w')

vcf = pysam.VariantFile(currentDirectory + '/experiment_0.vcf')
samples = np.array(vcf.header.samples)
variants = list(vcf.fetch())


#Question 1:
print('\nQuestion 1', file=outfile)
print('Number of people in the experiment: {}'.format(len(samples)), file=outfile)
infected_num = 0
for sample in samples:
    if sample.startswith('case'):
        infected_num += 1
print('Number of infected people: {}'.format(infected_num), file=outfile)
print('Number of healthy people: {}'.format(len(samples)-infected_num), file=outfile)

# Question 2:
print('\n\nQuestion 2',file=outfile)
print ('Number of mutations in total:', len(variants), file=outfile)
chrVariants = {}
for v in variants:
    if chrVariants.get(v.chrom) is None:
        chrVariants[v.chrom] = 0
    chrVariants[v.chrom]+=1
for chrm in chrVariants.keys():
    print('chromosome {} | mutations num {}'.format(chrm, chrVariants[chrm]),file=outfile)

#Question 3:
from collections import Counter

def check_freq(variant):
  counts = Counter([sum(sample['GT']) for sample in variant.samples.values() if None not in sample['GT']])
  total = sum(counts.values())
  p = (2*counts[0]+counts[1])/(2*total)
  q = (2*counts[2]+counts[1])/(2*total)
  ref_homoZ = counts[0]/total
  alt_homoZ = counts[2]/total
  heteroZ = counts[1]/total
  return min(p,q), ref_homoZ, alt_homoZ, heteroZ

mafArr = []
refArr = []
altArr = []
heteroArr = []
for variant in variants:
    maf, ref, alt, het = check_freq(variant)
    mafArr.append(maf)
    refArr.append(ref)
    altArr.append(alt)
    heteroArr.append(het)

plt.figure(figsize=(12,10))
plt.title('Question 3')
plt.subplot(221)


plt.hist(mafArr, bins=50, edgecolor='black', linewidth=1.0)
plt.xlabel('MAF')

plt.subplot(222)
plt.hist(refArr, bins=50, edgecolor='black', linewidth=1.0)
plt.xlabel('REF Alelle Homozygot')

plt.subplot(223)
plt.hist(altArr, bins=50, edgecolor='black', linewidth=1.0)
plt.xlabel('ALT Alelle Homozygot')

plt.subplot(224)
plt.hist(heteroArr, bins=50, edgecolor='black', linewidth=1.0)
plt.xlabel('Heterozygot')

plt.savefig(currentDirectory + '/3_plots.png')
print('\n\nQuestion 3\n Plots saved in a file 3_plots.png', file=outfile)

# Question 4:
from scipy.stats import chisquare


outfile.close()

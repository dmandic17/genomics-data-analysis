import os
import pysam
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


# Intro code (all questions need this): 
currentDirectory = os.path.dirname(os.path.realpath(__file__))
print(currentDirectory)
outfile = open(currentDirectory + '/generated_files/outfile.txt', 'w')

vcf = pysam.VariantFile(currentDirectory + '/experiment_0.vcf')
samples = np.array(vcf.header.samples)
variants = list(vcf.fetch())

case_all = {x for x in samples if x.startswith('case')}
control_all = {x for x in samples if x.startswith('control')}


# Question 1:
print('\nQuestion 1', file=outfile)
print('Number of people in the experiment: {}'.format(len(samples)), file=outfile)
print('Number of infected people: {}'.format(len(case_all)), file=outfile)
print('Number of healthy people: {}'.format(len(control_all)), file=outfile)


# Question 2:
print('\n\nQuestion 2',file=outfile)
print ('Number of mutations in total:', len(variants), file=outfile)
chrVariants = {}
for v in variants:
    if chrVariants.get('chr '+ v.chrom) is None:
        chrVariants['chr '+ v.chrom] = 0
    chrVariants['chr '+ v.chrom]+=1
print(pd.DataFrame.from_dict(chrVariants, orient='index', columns=['Num. of mutations']), file=outfile)

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

plt.savefig(currentDirectory + '/generated_files/3_plots.png')
print('\n\nQuestion 3\nPlots saved in a file 3_plots.png', file=outfile)


# Question 4:
#4a
from scipy.stats import chisquare

def run_chi2(variant):
  case_refs = sum(1 for sample in case_all for gt in variant.samples[sample]['GT'] if gt==0)
  case_alts = sum(1 for sample in case_all for gt in variant.samples[sample]['GT'] if gt==1)
  control_refs = sum(1 for sample in control_all for gt in variant.samples[sample]['GT'] if gt==0)
  control_alts = sum(1 for sample in control_all for gt in variant.samples[sample]['GT'] if gt==1)

  total = case_alts + case_refs + control_alts + control_refs
  
  expected = [(case_refs + case_alts)*(case_refs + control_refs)/total,
              (control_refs + control_alts)*(case_refs + control_refs)/total,
              (case_refs + case_alts)*(case_alts + control_alts)/total,
              (control_refs + control_alts)*(case_alts + control_alts)/total]
  return chisquare([case_refs, control_refs, case_alts, control_alts], expected).pvalue

print('\n\nQuestion 4:\n4a -> Running chi-square...', file=outfile)
chi2 = [(variant.chrom, variant.pos, run_chi2(variant)) for variant in variants]
print ('Done.', file=outfile)

#4b
print('4b -> The correction is neccesary, because we have over one milion mutations, and that effects the limit for p-value.', file=outfile)
print('I\'m using Bonferroni correction so as the limit I\'m using 0.05/len(data) in code I used len(chi2) but it is the same length.', file=outfile)
print('Bonferroni correction is just changing the limit, but there are also other corrections in use.', file=outfile)

#4c
from itertools import groupby

results = {key: list(value) for key, value in groupby(chi2, key=lambda x: x[0])}
remember = 0

plt.figure(figsize=(13, 10))
plt.subplot(211)
ticks = []

for i, chrom in enumerate(results):
  pos = [x[1] for x in results[chrom]]
  pos = np.array(pos)
  ticks.append(remember + max(pos)/2)
  pos += remember
  remember = max(pos)
  p = [-np.log10(x[2]) for x in results[chrom]]
  plt.scatter(pos, p)
  plt.axhline(-np.log10(0.05), c='r')  
tickLabels = [str(i) for i in range(1,23)]
plt.title('Manhattan diagram before correction')
plt.xticks(ticks, tickLabels)

plt.subplot(212)
ticks = []
for i, chrom in enumerate(results):
  pos = [x[1] for x in results[chrom]]
  pos = np.array(pos)
  ticks.append(remember + max(pos)/2)
  pos += remember
  remember = max(pos)
  p = [-np.log10(x[2]) for x in results[chrom]]
  plt.scatter(pos, p)
  plt.axhline(-np.log10(0.05/len(chi2)), c='r')  #correction here
tickLabels = [str(i) for i in range(1,23)]
plt.title('Manhattan diagram after correction')
plt.xticks(ticks, tickLabels)
plt.savefig(currentDirectory + '/generated_files/4c_manhattan.png')
print('4c -> Plots saved in a file 4c_manhattan.png', file=outfile)

#4d
# Hardy-Weinberg equilibrium:
from collections import Counter

def check_hwe(variant):
  counts = Counter([sum(sample['GT']) for sample in variant.samples.values() if None not in sample['GT']])
  total = sum(counts.values())
  p = (2*counts[0]+counts[1])/(2*total)
  q = (2*counts[2]+counts[1])/(2*total)

  observed = [counts[0], counts[1], counts[2]]
  expected = [p*p*total, 2*p*q*total, q*q*total]

  return chisquare(observed, expected, ddof=1).pvalue

print('4d -> Running HWE...', file=outfile)
hwe_fail = [check_hwe(variant) < 1e-6 for variant in variants]

chi2_filtered = [chi for chi, hwe in zip(chi2, hwe_fail) if not hwe]
results_filtered = {key: list(value) for key, value in groupby(chi2_filtered, key=lambda x: x[0])}
print('Done', file=outfile)

# Manhattan again after HWE:
plt.figure(figsize=(13, 5))
remember = 0
ticks = []
for i, chrom in enumerate(results):
  pos = [x[1] for x in results_filtered[chrom]]
  pos = np.array(pos)
  ticks.append(remember + max(pos)/2)
  pos += remember
  remember = max(pos)
  p = [-np.log10(x[2]) for x in results_filtered[chrom]]
  plt.scatter(pos, p)
  plt.axhline(-np.log10(0.05/len(chi2)), c='r')
tickLabels = [str(i) for i in range(1,23)]
plt.xticks(ticks, tickLabels)

plt.title('Manhattan after HWE')
plt.savefig(currentDirectory + '/generated_files/4d_manhattan.png')
print('Plot after hwe saved in a file 4d_manhattan.png', file=outfile)


# Question 5:
print('\n\nQuestion 5:', file=outfile)
print([('chromosome: '+str(val[0]),'position: '+ str(val[1])) for val in chi2 if val[2]<0.05/len(chi2)], file=outfile)
print ('UCSC Genome Browser: Chromosome 1, position 16896217 constains gene called NBPF1. (GRCh37/hg19)', file=outfile)
print('NBPF stands for neuroblastoma breakpoint family, member 1. ', file=outfile)
print('NBPF1 is a Protein Coding gene. Diseases associated with NBPF1 include Astrocytoma and Neuroblastoma. (GeneCard)', file=outfile)
print('More info on colab: https://colab.research.google.com/drive/1XEzV1tGHgnz4kWGpy_nA8XaSndYKHc4W#scrollTo=PUKNRrJq8VsI', file=outfile)
print('\n\nFinished.', file=outfile)
outfile.close()

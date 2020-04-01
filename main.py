import pysam
import numpy as np

vcf = pysam.VariantFile('experiment_0.vcf')
samples = np.array(vcf.header.samples)

# Question 1:
print('\nQuestion 1')
print('Number of people in the experiment: ', len(samples))
infected_num = 0
for sample in samples:
    if sample.startswith('case'):
        infected_num += 1
print('Number of infected people: ', infected_num)
print('Number of healthy people: ', len(samples)-infected_num)

# Question 2:
print('\nQuestion 2')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

# Initialize argument parser
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="Name of the input file")
parser.add_argument("output_file", help="Name of the output file")
args = parser.parse_args()

# list of genuine chromosomes
chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']

# load your file
df = pd.read_csv(args.input_file, sep='\t', header=None)

# Count the initial number of entries
initial_count = len(df)

# filter the data for only the genuine chromosomes and junction tag > 1
df = df[df[0].isin(chromosomes) & (df[3] > 1)]

# Count the number of entries after filtering
filtered_count = len(df)

# save the filtered data to the output file
df.to_csv(args.output_file, sep='\t', header=False, index=False)

# calculate length of each eccDNA
df['Length'] = df[2] - df[1]

# filter length from 0 to 2000
df = df[(df['Length'] >= 0) & (df['Length'] <= 2000)]

# calculate abundance of eccDNAs from each chromosome
chromosome_abundance = df[0].value_counts()

# calculate length distribution
length_distribution = df['Length'].value_counts().sort_index()


# create line plot for length distribution
plt.figure(figsize=(10, 6))
length_distribution.plot(kind='line')
plt.title('EccDNA Length Distribution')
plt.xlabel('Length')
plt.ylabel('Frequency')
plt.xticks(np.arange(0, 2000 + 1, 200))  # Set the x-ticks to be in increments of 200
plt.xlim([0, 2000])  # Set the limits of x-axis from 0 to 2000
plt.savefig('length_distribution.png')

# create bar plot for chromosome abundance
plt.figure(figsize=(10, 6))
chromosome_abundance.plot(kind='bar')
plt.title('EccDNA Chromosome Abundance')
plt.xlabel('Chromosome')
plt.ylabel('Frequency')
plt.savefig('chromosome_abundance.png')

# create bar plot for initial and filtered counts
plt.figure(figsize=(10, 6))
plt.bar(['Initial', 'Filtered'], [initial_count, filtered_count])
plt.title('Number of Entries Before and After Filtering')
plt.xlabel('Stage')
plt.ylabel('Count')
plt.savefig('filtering_counts.png')














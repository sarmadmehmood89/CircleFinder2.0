import pandas as pd
import matplotlib.pyplot as plt

# Paths for input files
path_eccdna = 'path/to/your/eccDNA/file.txt'
path_gtf = 'path/to/your/GTF/file.gtf'
path_introns = 'path/to/your/introns/file.bed'

# Load your data
eccdna_df = pd.read_csv(path_eccdna, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Junction_Tags'])

# Load GTF file
gtf_df = pd.read_csv(path_gtf, sep='\t', comment='#', header=None)
gtf_df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

# Load introns file
introns_df = pd.read_csv(path_introns, sep='\t', header=None)
introns_df.columns = ['seqname', 'start', 'end']
introns_df['feature'] = 'intron'

# Combine GTF and introns dataframes
all_features_df = pd.concat([gtf_df, introns_df], ignore_index=True)

# Filter to relevant features
features = ['exon', 'UTR', 'intron']
all_features_df = all_features_df[all_features_df['feature'].isin(features)]

# Create a dictionary to hold counts for each feature
feature_counts = {}

# Create a dictionary to hold sizes for each feature
feature_sizes = {}

# Loop over each feature
for feature in features:
    # Filter to rows of this feature
    feature_df = all_features_df[all_features_df['feature'] == feature]

    # Compute the size of this feature by summing the lengths of all instances
    feature_sizes[feature] = (feature_df['end'] - feature_df['start']).sum()

    # Loop over each row in eccdna_df
    for _, row in eccdna_df.iterrows():
        chrom, start, end, _ = row

        # Check if this region overlaps with the feature
        overlaps = ((feature_df['seqname'] == chrom) & (feature_df['start'] <= end) & (feature_df['end'] >= start)).any()
        if overlaps:
            feature_counts[feature] = feature_counts.get(feature, 0) + 1

# Calculate the percentage
feature_percentage = {feature: (count / len(eccdna_df)) * 100 for feature, count in feature_counts.items()}

# Print the percentage
for feature, percentage in feature_percentage.items():
    print(f'The percentage of eccDNAs coming from {feature}s is {percentage:.2f}%')

# Calculate fold enrichment
genome_size = 3.2e9  # Human genome size
feature_fold_enrichment = {feature: (count / len(eccdna_df)) / (feature_sizes[feature] / genome_size) for feature, count in feature_counts.items()}

# Print the fold enrichment
for feature, fold_enrichment in feature_fold_enrichment.items():
    print(f'The fold enrichment of eccDNAs from {feature}s is {fold_enrichment:.2f}')

# Draw a bar plot for percentage
plt.figure(figsize=(10, 6))
plt.bar(feature_percentage.keys(), feature_percentage.values())
plt.xlabel('Feature')
plt.ylabel('Percentage')
plt.title('Percentage of eccDNAs from each feature')
plt.savefig('genomic_feature_percentage.png')

# Draw a bar plot for fold enrichment
plt.figure(figsize=(10, 6))
plt.bar(feature_fold_enrichment.keys(), feature_fold_enrichment.values())
plt.xlabel('Feature')
plt.ylabel('Fold Enrichment')
plt.title('Fold Enrichment of eccDNAs from each feature')
plt.savefig('genomic_feature_fold_enrichment_hct_complete.png'))






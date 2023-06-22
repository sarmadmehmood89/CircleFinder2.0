import pysam
import matplotlib.pyplot as plt
import sys

# Function to calculate microhomology
def calculate_microhomology(sequence, k):
    left_end = sequence[:k]  # k bases from the left end
    right_end = sequence[-k:]  # k bases from the right end

    # Compare the nucleotides at each position and count microhomology
    microhomology_count = sum(1 for a, b in zip(left_end, right_end) if a == b)

    return microhomology_count > 0

# Input files from command line arguments
genome_file = sys.argv[1] 
eccdna_file = sys.argv[2]

# Open the genome file
genome = pysam.FastaFile(genome_file)

# Open the eccDNA file
with open(eccdna_file, "r") as file:
    lines = file.readlines()

# Counters for eccDNA regions with and without microhomology
with_microhomology = 0
without_microhomology = 0

# Process each line in the eccDNA file
for line in lines:
    chr_name, start, end, _ = line.strip().split()
    start, end = int(start), int(end)

    # Fetch the sequence from the genome
    sequence = genome.fetch(chr_name, start, end)

    # Check for microhomology
    if calculate_microhomology(sequence, k=15):  # Here, k is set to 15
        with_microhomology += 1
    else:
        without_microhomology += 1

total_eccDNA = with_microhomology + without_microhomology

# Calculate percentages
with_microhomology_percentage = (with_microhomology / total_eccDNA) * 100
without_microhomology_percentage = (without_microhomology / total_eccDNA) * 100

# Plot the percentages
plt.bar(['With microhomology', 'Without microhomology'],
        [with_microhomology_percentage, without_microhomology_percentage], color='b')

plt.xlabel('Category')
plt.ylabel('Percentage')  # y-axis label is now 'Percentage'
plt.title('Microhomology Analysis')

# Save the plot to a file before showing it
plt.savefig('microhomology_plot.png')

plt.show()

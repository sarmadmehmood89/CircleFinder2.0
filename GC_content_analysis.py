import pysam
import matplotlib.pyplot as plt
import pandas as pd

def calculate_gc_content(sequence):
    """Calculates and returns the GC content of the sequence."""
    return (sequence.count('G') + sequence.count('C')) / len(sequence) * 100

def main():
    genome_file = "/path_to/hg38.fa"
    eccdna_file = "/path_to/eccDNA.txt"
    
    # Try to open the genome file
    try:
        genome = pysam.FastaFile(genome_file)
    except IOError:
        print(f"Could not read file: {genome_file}")
        return

    # Try to open the eccDNA file
    try:
        with open(eccdna_file, "r") as file:
            lines = file.readlines()
    except IOError:
        print(f"Could not read file: {eccdna_file}")
        return

    gc_content_values = []

    for line in lines:
        chr_name, start, end, _ = line.strip().split()
        start, end = int(start), int(end)
        
        # Fetch the sequence from the genome
        sequence = genome.fetch(chr_name, start, end)
        
        # Calculate GC content and add it to the list
        gc_content = calculate_gc_content(sequence)
        gc_content_values.append(gc_content)

    # Convert the list into pandas DataFrame
    df = pd.DataFrame(gc_content_values, columns=['GC content'])

    # Plot a histogram
    df['GC content'].plot(kind='hist', edgecolor='black')
    plt.xlabel('GC Content (%)')
    plt.title('Histogram of GC Content across eccDNA regions')
    plt.savefig('gc_content_histogram.png')
    plt.show()

    # Plot a density plot
    df['GC content'].plot(kind='density')
    plt.xlabel('GC Content (%)')
    plt.title('Density Plot of GC Content across eccDNA regions')
    plt.savefig('gc_content_density.png')
    plt.show()

if __name__ == "__main__":
    main()

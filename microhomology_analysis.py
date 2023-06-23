import pysam
import matplotlib.pyplot as plt
import sys

def calculate_microhomology(seq1, seq2):
    """Checks for microhomology between two sequences"""
    min_length = min(len(seq1), len(seq2))
    for k in range(min(15, min_length), 1, -1):
        if seq1[:k] == seq2[:k]:
            return True
    return False

def get_sequences(genome, chr_name, start, end):
    """Fetches upstream and downstream sequences for a given genomic region"""
    upstream_start = genome.fetch(chr_name, start-15, start)
    upstream_end = genome.fetch(chr_name, end-15, end)
    downstream_start = genome.fetch(chr_name, start, start+15)
    downstream_end = genome.fetch(chr_name, end, end+15)
    return upstream_start, upstream_end, downstream_start, downstream_end

def main(genome_file, eccdna_file):
    genome = pysam.FastaFile(genome_file)

    with open(eccdna_file, "r") as file:
        lines = file.readlines()

    with_microhomology, without_microhomology = 0, 0

    for line in lines:
        chr_name, start, end, _ = line.strip().split()
        start, end = int(start), int(end)

        upstream_start, upstream_end, downstream_start, downstream_end = get_sequences(genome, chr_name, start, end)

        if calculate_microhomology(upstream_start, upstream_end) or calculate_microhomology(downstream_start, downstream_end):
            with_microhomology += 1
        else:
            without_microhomology += 1

    total_eccDNA = with_microhomology + without_microhomology

    with_microhomology_percentage = (with_microhomology / total_eccDNA) * 100
    without_microhomology_percentage = (without_microhomology / total_eccDNA) * 100

    plt.bar(['With microhomology', 'Without microhomology'],
            [with_microhomology_percentage, without_microhomology_percentage], color='b')

    plt.xlabel('Category')
    plt.ylabel('Percentage')  
    plt.title('Microhomology Analysis')
    plt.savefig('microhomology_plot.png')
    plt.show()

if __name__ == "__main__":
    genome_file = sys.argv[1]
    eccdna_file = sys.argv[2]
    main(genome_file, eccdna_file)

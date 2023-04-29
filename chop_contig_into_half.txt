#!/usr/bin/perl -w
use strict;

# Check if the number of input arguments is correct
die "Usage: perl $0 contig.fna \n" if (@ARGV != 1);

# Define the output directory
my $dir = "./chopContig";
system("mkdir chopContig");

# Load the necessary library
use Bio::SeqIO;

# Get input file name from command line arguments
my $input = shift @ARGV;

# Initialize Bio::SeqIO object to read the input file
my $in = new Bio::SeqIO(-format => 'fasta', -file => "$input");

# Initialize variables
my $count = 0;
my $len = 0;
my $middle = 0;
my $left = '';
my $right = '';
my $name = '';

# Iterate through each sequence in the input file
while (my $seq = $in->next_seq) {
    $count++;
    $name = $seq->id;
    $len = $seq->length;
    $middle = int($len / 2);
    $left = $seq->subseq(1, $middle);
    $right = $seq->subseq(($middle + 1), $len);

    # Open output file for writing and handle errors
    open(OUT, ">$dir/$name.fa") || die "can't open file: $!\n";
    print OUT ">$name.left.sequence\n$left\n>$name.right.sequence\n$right\n";
    close(OUT); # Close the output file
}

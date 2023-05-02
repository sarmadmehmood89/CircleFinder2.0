#!/usr/bin/perl

# This script processes FASTA files in a specified directory and outputs a combined FASTA file
# specifically for circular contig sequences.

# Check for correct number of command-line arguments
die "Usage: $0 final_fasta(Dir) > circularcontig.fna" unless (@ARGV == 1);

# Construct the path to the input directory
my $dir = "./"."$ARGV[0]";

# Open the input directory
opendir(DIR, $dir) || die "Can't open directory $dir\n";

# Read and process the filenames from the input directory
@array = ();
@array = readdir(DIR);

my $name = '';
foreach my $file (@array){
    # Skip if the filename does not end with ".final.fasta"
    next unless ($file =~ /^\S+.final.fasta$/);

    # Capture the part of the filename before ".fasta" and assign it to the $name variable
    if ($file =~ /^(\S+).fasta$/){
        $name = $1;
    }

    # Open and process each FASTA file
    open (FILE, "$dir/$file")||die "can't open file:$!\n";
    while(<FILE>){
        chomp;

        # Skip empty lines
        if (/^$/){next;}

        # Print sequence identifier with the $name variable
        if(/^>/){
            print ">$name\n";
        }
        # Print sequence data
        elsif(/^\S+$/){
            print "$_\n";
        }
    }

    # Close the file
    close FILE;
}

# Close the directory
closedir (DIR);

#!/usr/bin/perl
use strict;
use File::Basename;

# Check if the number of input arguments is correct
die "Usage: perl $0 <directory_name>\n" if (@ARGV != 1);

# Define the input directory from command line arguments
my $dir = "./$ARGV[0]";
opendir(DIR, $dir) || die "Can't open directory $dir\n";

my @array = ();
@array = readdir(DIR);

# Iterate through each file in the input directory
foreach my $file (@array) {
    next unless ($file =~ /^\S+.fa$/);

    my $name = '';
    if ($file =~ /^(\S+).fa$/) {
        $name = $1;
    }

    # Convert fasta files to AMOS message format
    system("/path/to/amos/bin/toAmos -s $dir/$name.fa -o $dir/$name.afg");

    # Create a bank from the AMOS message file
    system("/path/to/amos/bin/bank-transact -c -z -b $dir/$name.bnk -m $dir/$name.afg");

    # Extract reference sequences and query sequences from the bank
    system("/path/to/amos/bin/dumpreads $dir/$name.bnk -M 1 > $dir/$name.ref.seq");
    system("/path/to/amos/bin/dumpreads $dir/$name.bnk -m 1 > $dir/$name.qry.seq");

    # Run nucmer to align reference and query sequences
    system("/path/to/MUMmer3.23/nucmer -c 10 $dir/$name.ref.seq $dir/$name.qry.seq -p $name");

    # Move the output delta file to the specified directory
    system("mv $name.delta $dir/");

    # Extract the alignment coordinates from the delta file
    system("/path/to/MUMmer3.23/show-coords -H -c -l -o -r -I 94 $dir/$name.delta | /path/to/amos/bin/nucmerAnnotate | egrep 'BEGIN|END|CONTAIN|IDENTITY' > $dir/$name.coords");

    # Convert the nucmer coordinates to AMOS overlaps
    system("/path/to/amos/bin/nucmer2ovl -ignore 10 -tab $dir/$name.coords | /path/to/amos/bin/sort2 > $dir/$name.ovl");

    
    # Check for the existence of the overlap file
    my $ovl_file = "$dir/$name.ovl";
    if (-e $ovl_file) {
        open(my $fh_ovl, '<', $ovl_file) or die "Could not open file '$ovl_file': $!";
        my $max_value = 0;

        # Find the maximum value in the 6th column of the overlap file
        while (my $line = <$fh_ovl>) {
            my @columns = split /\t/, $line;
            my $value = $columns[5];
            if ($value > $max_value) {
                $max_value = $value;
            }
        }
        close $fh_ovl;

        # If the maximum value is greater than 0, create a final.fasta file with the modified sequence
        if ($max_value > 0) {
            open(my $fh_in, '<', "$dir/$name.fa") or die "Could not open file '$dir/$name.fa': $!";
            open(my $fh_out, '>', "$dir/$name.final.fasta") or die "Could not create file '$dir/$name.final.fasta': $!";

            my $seq = "";
            my $lab = "";
            my $merged_seq = "";
            my $seq_count = 0;

            # Read the input file and merge the sequences
            while (<$fh_in>) {
                if (/^>/) {
                    if ($seq ne "") {
                        $merged_seq .= $seq;
                        $seq_count++;
                        $seq = "";
                    }
                    $lab = $_ if $seq_count == 0;
                } else {
                    chomp;
                    $seq .= $_;
                }
            }
            $merged_seq .= $seq;

            # Remove the maximum overlap value from the end of the merged sequence
            $merged_seq = substr($merged_seq, 0, length($merged_seq) - $max_value);

            # Write the modified sequence to the output file
            print $fh_out $lab;
            print $fh_out $merged_seq, "\n";

            close($fh_in);
            close($fh_out);
        }
    }
}


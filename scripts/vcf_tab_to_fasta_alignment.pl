#!/usr/bin/env perl

# Program to convert output of VCFtools' vcf-to-tab
# to FASTA alignment.

# Sample input file
#	$ head input.vcf.tab
#	chr10	94051	C	./	./	./	./	./	T/T
#	chr10	94056	T	./	./	./	./	./	C/C
#	chr10	94180	G	./	A/A	./	./	./	./

use strict;
use warnings;
use Getopt::Long;

my $exclude_het = 0;	# Should we exclude heterozygous SNPs?
my $exclude_missing = 0;	# Should we exclude missing?
my $output_ref  = 0;	# Should we output the reference calls?
my $debug;
my $input_tab;

my $usage = "usage: $0 [--exclude_het] [--output_ref] -i input.tab";

my $result = GetOptions (	"exclude_het"	=> \$exclude_het,
							"exclude_missing"	=> \$exclude_missing,
							"output_ref"	=> \$output_ref,
							"debug"	=> \$debug,
							"i=s"			=> \$input_tab
						) or die "Incorrect usage. $usage\n";

my $starting_col = 3;
if ($output_ref) {
	print STDERR "Including reference sequence. Remove --output_ref flag to exclude.\n";
	$starting_col = 2;
}

my %iupac = (
			'G/G' => 'G',
			'C/C' => 'C',
			'T/T' => 'T',
			'A/A' => 'A',

			'G/T' => 'K',
			'T/G' => 'K',
			'A/C' => 'M',
			'C/A' => 'M',
			'C/G' => 'S',
			'G/C' => 'S',
			'A/G' => 'R',
			'G/A' => 'R',
			'A/T' => 'W',
			'T/A' => 'W',
			'C/T' => 'Y',
			'T/C' => 'Y',

			'./.' => '.',
		);

open (TAB, "<$input_tab")
	or die "ERROR: Could not open input file $input_tab.\n";

my $header = <TAB>;

my @col_names = split /\t/, $header;

# Make temporary file with just lines we're going to use
my $temp_tab = $input_tab . "_clean";
open (TEMP, ">$temp_tab")
	or die "ERROR: Could not open temp file $temp_tab.\n";

# Get number of columns
my $num_cols = scalar @col_names;
print STDERR "Number of columns:\t$num_cols\n";
my $COUTLS=0; # count total lines
my $COUALS=0; # count accepted lines

LINE: foreach my $line (<TAB>) {
	$COUTLS++;

	my @data = split /\t/, $line;
	
	# Skip if this is indel (Length of @data will be less than $num_cols)
	if ((scalar @data) < $num_cols) {
		print STDERR "Skipping indel.\n";
		next LINE;
	}
	
	# Skip if any basepairs are actually 2 or more together
	for (my $i = $starting_col; $i < $num_cols; $i++) {
		
		my $bp = $data[$i]; 
		chomp $bp;
		if ($bp =~ /\w{2,}/) { # two or more contiguous alphanumerics ...
			print STDERR "Skipping multi-basepair insertion.\n";
			next LINE;
		}
	}

	if ($exclude_het) {
		# Exclude heterozygotes. Keep only fixed SNPs
		for (my $i = $starting_col; $i < $num_cols; $i++) {
			
			my $bp = $data[$i]; 
			chomp $bp;
			if ($bp =~ /(\w)\/(\w)/) {
				if ($1 ne $2) {
					# commenting these out .. it works well. i.e. if there is a heterozygote, the entire snp position partakes no more.
					if($debug) {
						print STDERR "skipping het at $COUTLS: $1 vs. $2 @ colnum $i (line is \"$line\") | ";
					}
					next LINE;
				}
			}
		}
	}
	if ($exclude_missing) {
		for (my $i = $starting_col; $i < $num_cols; $i++) {
			
			my $bp = $data[$i]; 
			chomp $bp;
			if ($bp =~ /\.\/\.*/) {
				if($debug) {
					print STDERR "missing at $COUTLS: $1 vs. $2 @ colnum $i (line is \"$line\") | ";
				}
				next LINE;
			}
		}
	}
	
	# Otherwise write line to pure temporary file
	print TEMP $line;
	$COUALS++;
}
	
close TAB;
close TEMP;
print STDERR "$COUALS out of $COUTLS accepted: clean tab file at $temp_tab.\n";

# Now convert cleaned tabular file to FASTA alignment
my $seq = "";
my $b_ref = 0;

# we're going to go down the columns
for (my $i = $starting_col; $i < $num_cols; $i++) {

	my $ind = $col_names[$i];
	chomp $ind;
	
	print ">" . $ind . "\n";
	
	open (TEMP, "<$temp_tab")
		or die "ERROR: Could not open temp file $temp_tab.\n";

	# Count number of bp printed so far in this line
	my $count = 0;
	my $tcount = 0; # total counts
	
	foreach my $line (<TEMP>) {
	
		my @data = split /\t/, $line;
		
		my $nuc = $data[$i]; # $i is the column we're interested in
		chomp $nuc;
	
		if (! $b_ref){
			$seq .= $data[2];
		}
		# Infer and print basepair. There are a few possibilities 
		
		# If we're reference, just print basepair
		if ($i == 2) {
			print $nuc;
			$count++;
		
		# Haploid
		} elsif ($nuc =~ /(\w)\/$/) {
			print $1;
			$count++;
				
		# Missing data? Shouldn't have happened here.
		} elsif(($nuc eq './') || ($nuc eq './.')) {
			die "ERROR: Mssing data should not have reached this stage";
			$count++;
		
		# Data
		} elsif ($nuc =~ /(\w)\/(\w)/) {
			my $first = $1;
			my $second = $2;
			
			# Homozygote
			if ($first eq $second) {
				print $first;
				$count++;
			
			# Heterozygote: should have been got previously, at firststage
			} else {
				my $gt = $first . '/' . $second;
				print STDERR "2ndstage het: $gt\n";
				if ( !exists($iupac{$gt}) ) { die "ERROR: BP is $nuc\n"; }
                $gt = $iupac{$gt};
				print $gt;
				$count++;
			}
		} else {
			print STDERR "Not counted: linenum $tcount: $line\n"
		}
			
		if(($count !=0) && ($count%100 == 0)) {
			print "\n";
			# $count = 0;
		}
		$tcount++;
	}
	
	$b_ref = 1;
	close TEMP;
	print STDERR "$ind: $count accepted counts out of $tcount\n";
	
	print "\n";
}

# we're going to print the ref last
print ">ref\n";
$seq =~ s/(.{0,70}(.|$))/$1\n/g;
print $seq;

exit;

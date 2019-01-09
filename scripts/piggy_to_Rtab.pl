#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Pod::Usage;

# Convert piggy file to Rtab.
# assumes multi-threshold files are sorted by family/threshold (for excluding families/alleles on their initial frequency).

=head1  SYNOPSIS

 piggy_to_Rtab.pl -i /path/to/piggy.*.tab -o /path/to/output_file

 Input-Output:
 -i|--input		input csv file [required]
 -o|--output		output treeWAS input file [required]
 
 Gene frequency:
 -l|--low		min allele frequency to include in output 
			[default: 0.00]
 -h|--high		max allele frequency to include in output 
			[default: 1.00]
 
 Filtering options:
 -s|--samples		tab delimited list of samples to include in output 
 			[default: off]			
 -d|--dosage		upper threshold of dosage to exclude from alignment 
 			[default: 0]			
 -a|--all		include family cluster when processing files
			containing multiple alleles
			[default: exclude] 
			
 General options:
 -h|--help		usage information
 
=cut

# command line options
my $input = ''; 
my $output_file = '';

my $l_threshold = 0.0;
my $h_threshold = 1.00;

my $dosage_threshold = 0;

my $include_family = 0;
my $list = "";

my $help = 0;

GetOptions(

	'help|?' 	=> \$help,
	'input=s' 	=> \$input,
	'output=s'	=> \$output_file,
	'low=f' => \$l_threshold,
	'high=f' => \$h_threshold,
	'dosage=f' => \$dosage_threshold,
	'all' => \$include_family,
	'samples=s' => \$list,
	
) or pod2usage(1);
pod2usage(1) if $help;
pod2usage( {-message => q{input pirate file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input eq ''; 
pod2usage( {-message => q{output file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $output_file eq ''; 

# [optional] open list file
my %list  = (); 
my $no_samples_list = 0;
if ($list ne ''){
	open LIST, $list or die " - ERROR: could not open list ($list)\n";
	while (<LIST>){
	
		my $line = $_;
		chomp $line;
		
		my @vars = split(/\t/, $line, -1);
		
		$list {$vars[0]} = 1;
		
	}close LIST;
	
	# feedback
	$no_samples_list = keys(%list); 
	print " - $no_samples_list samples to include from list file.\n";
}

# parse input file.
my @headers = ();
my @samples = ();
my $no_samples = "";
my $idx = 14;

my %prop_thresh = ();
my %out_hash = ();
my %exclude = ();
my @include = ();
my %all = ();

my $prop_remove = 0;
my $dose_remove = 0;

my $gene_count = 0;
open INPUT, $input or die "Input file did not open.\n";
while(<INPUT>){
	
	my $line = $_;
	chomp $line;
	
	$line = substr $line, 1; # remove first "
	$line =~ s/"$//;
	
	my @line = split(/\"\,\"/, $line, -1);

	# get genome names
	if( /^\"Gene/ ){
		
		# check for samples in list 
		@headers = @line;	
		if ($list ne ''){
			for my $t ( $idx..$#headers ){ 
				if ($list{$headers[$t]}){
					$list{$headers[$t]} = 2;
					push(@include, $t);
				}
			}
		
		}else{
			
			for my $t ( $idx..$#headers ){ push(@include, $t)}; 	
			
		
		}
 
		# list missing samples
		my @missing_samples = ();
		for ( keys %list ){
			push(@missing_samples, $_) if $list{$_} == 1;
		}
		print " - missing samples:\n" if scalar(@missing_samples) > 0;
		print sprintf("%s\n", join("\n", @missing_samples)) if scalar(@missing_samples) > 0;
		
		# feedback and store
		@samples = @headers[@include];
		$no_samples = @include; 
		if ($list ne "" ){
			print " - $no_samples/$no_samples_list samples found in input file\n";
		}else{
			print " - $no_samples samples in input file\n";
		}
		
	}else{
	
		++$gene_count;
	
		# sanity check 
		die " - ERROR: header not found in file" if scalar(@headers) == 0; 
			
		# variables		
		my $a_name = $line[0];
		my $dosage = $line[5];
		
		# store presence/absence
		my $a_count = 0;
		for my $i (@include){
			$out_hash{$a_name}{$headers[$i]} = 1 if $line[$i] ne "";
			++$a_count if $line[$i] ne "";
		}
		
		# calculate proportion (all genomes) presence
		my $prop = $a_count/$no_samples;
		
		# exclude variants not within threshold frequencies
		my $exclude = 0;
		if ( ($prop < $l_threshold) || ($prop > $h_threshold) ){
			$exclude = 1;
			$prop_remove++;
		}
		if ( ($dosage > $dosage_threshold) && ($dosage_threshold > 0) ){
			$exclude = 1;
			$dose_remove++;
		}
		
		if ( $exclude == 0 ){
			$prop_thresh{$a_name} = $prop;
		}
	}
		
}close INPUT;

# feedback
my $no_included = keys %prop_thresh;
print " - IGR freq. thresholds: $l_threshold - $h_threshold\n";
print " - dosage threshold: < $dosage_threshold\n" if $dosage_threshold != 0;
print " - $prop_remove IGRs were not between frequency thresholds\n";
print " - $dose_remove IGRs were greater than dosage threshold\n";;
print " - $no_included of $gene_count IGRs were included\n";

# print to file
open OUT, ">$output_file" or die " - ERROR: output file ($output_file) would not open for writing\n";

# identify genes to include
my @included = sort(keys(%prop_thresh));

# headers
print OUT join("\t", "Gene", @samples ), "\n"; 

# per sample
for my $a (@included){
	
	# print binary present/absent
	my @outline = ($a);
	
	for my $sample(@samples){
	
		if ($out_hash{ $a }{ $sample }){
			push(@outline, "1" );
		}else{
			push(@outline, "0" );
		}
		
	}
	
	# print 
	print OUT join("\t",@outline), "\n";
	
}close OUT;

exit

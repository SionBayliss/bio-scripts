#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Pod::Usage;

# create piggy sequence alignment and gff

=head1  SYNOPSIS

 piggy_align.pl -i /path/to/piggy_folder -o /path/to/output_dir

 Input-Output:
 -i|--input		input csv file [required]
 -o|--output		output directory [required]
 
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
			
 General options:
 -p|--prefix	prefix for output file [default: IGR_alignment]
 -e|--exclude-off include a sequence from isolates which contain multiple sequences
 			for a single IGR (default: exclude and replace with dashes)
 -h|--help		usage information
 
=cut

# command line options
my $input = ''; 
my $output = '';

my $l_threshold = 0.0;
my $h_threshold = 1.00;

my $dosage_threshold = 0;

my $include_family = 0;
my $list = "";
my $exclude = 1; 
my $prefix = "IGR_alignment";

my $help = 0;

GetOptions(

	'help|?' 	=> \$help,
	'input=s' 	=> \$input,
	'output=s'	=> \$output,
	'low=f' => \$l_threshold,
	'high=f' => \$h_threshold,
	'dosage=f' => \$dosage_threshold,
	'samples=s' => \$list,
	'exclude-off' => $exclude,
	'prefix=s'	=> \$prefix,
	
) or pod2usage(1);
pod2usage(1) if $help;
pod2usage( {-message => q{input piggy dir is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input eq ''; 
pod2usage( {-message => q{output file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $output eq ''; 

# check output directory exists or make it
unless (-d $output ){
	mkdir($output) or die " - ERROR: could not make $output directory.\n";
}

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

my @include = ();
my @cluster_array = ();

my $prop_remove = 0;
my $dose_remove = 0;

my $gene_count = 0;
open INPUT, "$input/IGR_presence_absence.csv" or die "$input/IGR_presence_absence.csv file did not open.\n";
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
			push(@cluster_array, $a_name);
		}
	}
		
}close INPUT;

# feedback
my $no_included = @cluster_array;
print " - IGR freq. thresholds: $l_threshold - $h_threshold\n";
print " - dosage threshold: < $dosage_threshold\n" if $dosage_threshold != 0;
print " - $prop_remove IGRs were not between frequency thresholds\n";
print " - $dose_remove IGRs were greater than dosage threshold\n";;
print " - $no_included of $gene_count IGRs were included\n";

# create temp fasta files
mkdir("$output/isolate_core_IGR_tmp/"); 
for my $s (@samples){
	open OUTPUT_TMP, ">$output/isolate_core_IGR_tmp/$s.fasta" or die "Cannot open output file: $output/isolate_core_IGR_tmp/$s.fasta\n";
	print OUTPUT_TMP ">$s\n";
}

my @lengths = ();
my $alignment_length = 0;
foreach my $cluster(@cluster_array){
	
	my %cluster_seq_hash=();
	my $len="";
	
	my $cluster_id="";
	my @cluster_id_array=();
	my $isolate="";
	
	my %seq_count = ();
	
	open INPUT, "$input/cluster_intergenic_alignment_files/${cluster}_aligned.fasta" or die "Input file doesn't exist: $input/cluster_intergenic_alignment_files/${cluster}_aligned.fasta\n";
	while(my $line=<INPUT>){
	
		chomp $line;

		if($line =~ /^>(.+)/){
			$cluster_id = $1;
			@cluster_id_array = split(/_\+_\+_/, $cluster_id);
			$isolate = $cluster_id_array[0];
		}else{
		
			my $seq = uc($line);
			$seq =~ s/-/N/g;
			
			# increment sequence count per isolate
			$seq_count{$isolate}++;
			
			# store sequence
			$cluster_seq_hash{$isolate}=$seq;
			
			my $len_new=length($seq);
			
			# sanity check
			if ($len eq ""){
				$len = $len_new;
			}elsif($len != $len_new){
				print " - warning: length of $cluster_id in $isolate do not match reference length\n";	
			}

		}
	}
	close INPUT;
	
	foreach $isolate(@samples){
		open OUTPUT_TMP, ">>$output/isolate_core_IGR_tmp/$isolate.fasta" or die "Cannot open output file: $output/isolate_core_IGR_tmp/$isolate.fasta\n";
		
		# no sequence (Ns)
		if(!$cluster_seq_hash{$isolate}){
			for(my $i=0; $i<$len; $i++){
				print OUTPUT_TMP "N";
			}
		}
		# [optional] exclude sequence due to copy no. > 1 (-s)
		elsif( ($seq_count{$isolate} > 1) && ($exclude == 1) ){ 
			for(my $i=0; $i<$len; $i++){
				print OUTPUT_TMP "-";
			}
		}
		# otherwise print
		else{
			print OUTPUT_TMP "$cluster_seq_hash{$isolate}";
		}
	}
	
	# add length to array for gff
	push(@lengths, $len);
	$alignment_length += $len;
}
close OUTPUT_TMP;

foreach my $isolate(@samples){
	open OUTPUT_TMP, ">>$output/isolate_core_IGR_tmp/$isolate.fasta" or die "Cannot open output file: $output/isolate_core_IGR_tmp/$isolate.fasta\n";
	print OUTPUT_TMP "\n";
}
close OUTPUT_TMP;

# make gff
open GFF, ">$output/$prefix.$l_threshold-$h_threshold.gff" or die "Cannot open output file: $output/$prefix.$l_threshold-$h_threshold.gff\n";

# headers 
print GFF "##gff-version 3\n##sequence-region piggy 1 $alignment_length\n";

my $rolling = 0;
foreach my $i (0..$#cluster_array){

	my $cluster_id = $cluster_array[$i];
	my $cluster_len = $lengths[$i];
	
	# current start
	my $start = $rolling+1;
	my $end = $start+($cluster_len-1);
	
	# recalc rolling
	$rolling = $end;
	
	# make gff line for cluster
	my $gff_line = sprintf( "Piggy\tNA\tIGR\t%s\t%s\t\.\t\+\t0\tID=%s" , $start , $end , $cluster_id );
	print GFF "$gff_line\n";
}
close GFF;
print " - IGR GFF3 created.\n";

# open output alignmnet.
open OUTPUT, ">$output/$prefix.$l_threshold-$h_threshold.fasta" or die "Cannot open output file: $output/$prefix.$l_threshold-$h_threshold.fasta\n";
foreach my $isolate(@samples){
	open INPUT, "$output/isolate_core_IGR_tmp/$isolate.fasta" or die "Input file doesn't exist: $output/isolate_core_IGR_tmp/$isolate.fasta\n";
	while(my $line=<INPUT>){
		print OUTPUT "$line";
	}
}
close OUTPUT;
print " - IGR alignment created ($alignment_length bp).\n";

# remove temp
foreach my $isolate (@samples){
	unlink("$output/isolate_core_IGR_tmp/$isolate.fasta");
}
rmdir("$output/isolate_core_IGR_tmp/");



#!/usr/bin/env perl

use strict;
use warnings qw(all);

use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;
use File::Basename;
use Cwd 'abs_path';

# prepare data for pyseer

=head1  SYNOPSIS

	prepare_data_for_pyseer 
	
 Input variants:
 --scripts	path to PIRATE root installation [required] 
 --meta	metadata file [required] 
 --output 	output directory [required] 
 --pirate 	input PIRATE directory [optional]
 --vcf	input vcf file [optional]
 [required: either --vcf or --pirate]
 
 PIRATE options:
 --core	frequency threshold to be considered core [default: 0.95]
 --dosage	dosage threshold for inclusion [default: 1.2]
 --family-low	lower bounds threshold for including accessory genes [default: 0.05]
 
 VCF options:
 --prefix-vcf	set output filename prefix for vcf files [default: vcf]
 
 -h|--help 	usage information
  
=cut


# make command line 
my $commandline = $0 . " ". (join " ", @ARGV);

# command line options
my $vcf = "";
my $pir = "";
my $path = "";
my $phenotype = '';
my $output = '';
my $af = "";

my $core = "0.95"; 
my $dosage = "1.2";
my $family_low = "0.5";

my $prefix_vcf = "vcf"; 

my $help = 0;


pod2usage(1) if scalar(@ARGV) == 0;
GetOptions(

	'vcf=s'	=> \$vcf,
	'pirate=s'	=> \$pir,
	
	'metadata=s' => \$phenotype,
	'output=s' => \$output,
	'scripts=s' => \$path,
	
	'af=f' => \$af,
	
	'core=f' => \$core,
	'dosage=f' => \$dosage,
	'family-low=f' => \$family_low,
	
	'prefix-vcf=s' => \$prefix_vcf,
	
	'help|?' 	=> \$help,

) or die pod2usage(1);
pod2usage(1) if $help;


# check options
die " - ERROR: must specify phenotype file with --metadata\n" if $phenotype eq "";
die " - ERROR: must specify output directory with --output\n" if $output eq "";
die " - ERROR: must specify PIRATE or VCF file with --vcf/--pirate\n" if ( ($vcf eq "") && ($pir eq "") );
die " - ERROR: must specify PIRATE script path --path\n" if ( ($path eq "") && ($pir ne "") );
die " - ERROR: must specify allelic frequency with --af\n" if $af eq "";

# find path of bio-scripts 
my $script_path = abs_path(dirname($0));
my $pirate_path = abs_path($pir);

# set thresholds
die " - ERROR: allelic frequency must be between 0-1.\n" if ( ($af>1) && ($af<1) );
my $af_low = $af;
my $af_high = 1-$af;

# make log file
my $log = "$output/log.$af.txt";
system("echo $commandline > $log");

# create sample list
my $s_list = "$output/sample_list.$af.txt";
open META, "$phenotype" or die " - ERROR: could not open phenotype file ($phenotype)\n"; 
open SLIST, ">$s_list" or die " - ERROR: could not write to $s_list\n"; 
my $s_count = 0;
my $count = 0;
while(<META>){
	++$count;
	
	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line, -1);
	
	if ($count > 1){
		print SLIST "$vars[0]\n";
		++$s_count;
	}
	
}close SLIST;
print "\n - $s_count samples in metadata file\n";
system("echo - $s_count samples in metadata file >> $log");

# make log subroutine
sub print_log{
	my $val = shift(@_);
	print "$val";
	open LOG, ">>$log" or die " - ERROR: could not open log file\n";
	print LOG $val;
	close LOG;
} 

# process PIRATE files
if ($pir ne ""){

	print_log "\n -------------------------- \n\n - Processing PIRATE:\n";

	# check for scripts/files in path
	my $conv_path = sprintf("%s/tools/convert_format/", $path); 
	die " - ERROR: could not find $conv_path/paralogs_to_Rtab.pl" unless ( -f  "$conv_path/paralogs_to_Rtab.pl");
	die " - ERROR: could not find $conv_path/PIRATE_to_Rtab.pl" unless ( -f  "$conv_path/PIRATE_to_Rtab.pl");
	die " - ERROR: could not find $pirate_path/PIRATE.unique_alleles.tsv" unless ( -f  "$pirate_path/PIRATE.unique_alleles.tsv");
	die " - ERROR: could not find $pirate_path/PIRATE.gene_families.tsv" unless ( -f  "$pirate_path/PIRATE.gene_families.tsv");
	
	# core_alleles
	print_log "\n - core alleles:\n\n";
	system("$conv_path/PIRATE_to_Rtab.pl -i $pirate_path/PIRATE.unique_alleles.tsv -o $output/core_alleles.$af.tsv --low $af_low --high $af_high -fl $core -d $dosage -s $s_list | tee -a $log");
	die " - ERROR: core alleles failed\n" if $?;
	
	# accessory_pe
	print_log "\n - accessory presence_absence:\n\n";
	system("$conv_path/PIRATE_to_Rtab.pl -i $pirate_path/PIRATE.gene_families.tsv -o $output/accessory_pa.$af.tsv --low $af_low --high $af_high -fl $family_low -fh $core -d $dosage -s $s_list | tee -a $log");
	die " - ERROR: accessory pa failed\n" if $?;
	
	# accessory_alleles
	print_log "\n - accessory alleles:\n\n";
	system("$conv_path/PIRATE_to_Rtab.pl -i $pirate_path/PIRATE.unique_alleles.tsv -o $output/accessory_alleles.$af.tsv --low $af_low --high $af_high -fl $family_low -fh $core -d $dosage -s $s_list | tee -a $log");
	die " - ERROR: accessory alleles failed\n" if $?;
	
	# duplications
	print_log "\n - duplications:\n\n";
	system("$conv_path/paralogs_to_Rtab.pl -i $pirate_path/PIRATE.unique_alleles.tsv -o $output/duplications.$af.tsv --low $af_low --high $af_high -b -t d -s $s_list | tee -a $log");
	die " - ERROR: duplications failed\n" if $?;
	
	# ff
	print_log "\n - fission fusion:\n\n";
	system("$conv_path/paralogs_to_Rtab.pl -i $pirate_path/PIRATE.unique_alleles.tsv -o $output/fissions.$af.tsv --low $af_low --high $af_high -b -t ff -s $s_list | tee -a $log");
	die " - ERROR: fission fusion failed\n" if $?;

}

# process vcf file
if ($vcf ne ""){

	print_log "\n -------------------------- \n\n - Processing VCF:\n";
	
	# check for scripts/files in path
	die " - ERROR: could not find $vcf" unless ( -f  "$vcf");
	
	# expand vcf  - filter low freq variants (dataset dependent)
	print_log "\n - expand vcf:\n";
	system("$script_path/expand_vcf -i $vcf -o $output/$prefix_vcf.expanded.$af.vcf -fl $af_low -fh $af_high -l $s_list | tee -a $log"); 
	die " - ERROR: expand vcf failed\n" if $?;
	
	# convert vcf to fasta (filter sites)
	print_log "\n - vcf to fasta:\n";
	system("$script_path/vcf_to_fasta -i $vcf -o $output/$prefix_vcf.filtered.$af.fasta -fl $af_low -fh $af_high -l $s_list | tee -a $log"); 
	die " - ERROR: vcf to fasta failed\n" if $?;

	# convert fasta to sim-mat
	print_log "\n - fasta to similarity matrix:\n\n";
	system("$script_path/snp-sim $output/$prefix_vcf.filtered.$af.fasta -b > $output/$prefix_vcf.simmat 2>>$log"); 
	die " - ERROR: fasta to similarity matrix failed\n" if $?;
	
}

# remove redundant lines from log file
system("grep -v 'variant sites processed' < $log > $log.temp");
system("mv $log.temp $log");

# make combined output
print_log "\n - combining outputs...\n\n";
system("$script_path/combine_variants_for_pyseer -v $output/$prefix_vcf.expanded.$af.vcf -p $output/core_alleles.$af.tsv -p $output/accessory_pa.$af.tsv -p $output/accessory_alleles.$af.tsv -p $output/duplications.$af.tsv -p $output/fissions.$af.tsv --suffix-pres \"c,a,aa,d,ff\" --samples $s_list -o $output/combined_variants.$af.tab | tee -a $log"); 
die " - ERROR: combining outputs failed\n" if $?;

# feedback
print_log " - complete\n\n";



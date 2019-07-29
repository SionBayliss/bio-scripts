#!/usr/bin/env perl

# Add what a gene encodes resistance to from abricate output.

use File::Basename;
use Cwd 'abs_path';
my $path = dirname(abs_path($0));

# Inputs
$table=$ARGV[0];

# Parse notes file in abricate dir.
print $path;
open NOTES, "$path/resfinder_notes.txt" or die $1;
while(<NOTES>){
	unless(/^#/){
		@entry=split(/:/,$_);
		$r{$entry[0]}=$entry[1];
		
		if($_=~/Alternate/){
			@entry2=split(/;/,$entry[2]);
			$r{$entry2[0]}=$entry2[1];
		}
	}	
}

# Parse table and add resistance column.
open TABLE, $table or die $!;
while(<TABLE>){
	$line=$_;
	chomp $line;
	
	# Change for plotting
	$line=~s/\.fsa||\.fasta||\.fa//g;
	$line=~s/#//g;
	
	if(/#FILE/){
		print "$line\tResistance\n"; 
	}else{
		@line=split(/\t/,$line);
		$gene=$line[4];
		$gene =~ s/\_\d+$//g;
		
		$gene = $1 if ($gene =~ /(mcr-\d)\.\d+/);
		
		unless(exists $r{$gene}){
			print "ERROR: $gene not found.\n";
		}else{
			$res=$r{$gene};
			$res=~s/\sresistance//g;
			print "$line\t", $res,"\n";
		}
	}
}close TABLE;

exit

#!/usr/bin/perl

use strict;
use warnings;

### SET GLOBAL PATHS -------------------------------------------------
my $dir_work = "/Data/Raunak/projects/MESO_peritoneal/data/cnv/";
my $dir_data = $dir_work."mes_wdpm/";

### DEFINE FILES ------------------------------------------------------
my $file_dat = $dir_data."meso_wdpm_seq_calls_matrix.txt";
my $file_prb = $dir_data."meso_wdpm_cnv_probe_medians_calls.tsv";
my $file_seg = $dir_data."meso_wdpm_cnv_seg_values_calls.tsv";

#print $file_dat, "\n";

### READ FILE ---------------------------------------------------------
open(FILE_DAT, $file_dat) or die "Error: Cannot open data file";
open(FILE_PRB, ">$file_prb") or die "Error: Cannot open output data file";
open(FILE_SEG, ">$file_seg") or die "Error: Cannot open output data file";

while(<FILE_DAT>){
	chomp;
	my $line = $_;
	my @fields = split("\t", $line);
	
	my $ary_size = @fields - 1;
	#my $ary_size = 80058;
	
	#print $ary_size, "\n";
	
	my $sample_name = $fields[0];
	
	#print $sample_name, "\n";

	print FILE_PRB $sample_name;
	print FILE_SEG $sample_name;
	
	### Probe median _odd_set_
	for(my $i = 1; $i <= $ary_size; $i = $i + 3){
		#print "\t", $fields[$i];
		print FILE_PRB "\t", $fields[$i];
	}
	print FILE_PRB "\n";
	
	### SEgment values _odd_set_
	for(my $i = 2; $i <= $ary_size; $i = $i + 3){
		print FILE_SEG "\t", $fields[$i];
	}
	print FILE_SEG "\n";
}
close FILE_DAT;
close FILE_PRB;
close FILE_SEG;

exit;


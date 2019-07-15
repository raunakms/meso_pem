#!/usr/bin/perl

use strict;
use warnings;

### SET GLOBAL PATHS -------------------------------------------------
my $dir_work = "/collinsgroup/Raunak/projects/MESO_peritoneal/data/cnv/";
my $dir_data = $dir_work."seq_call_refseq_genes_meso/";

### DEFINE FILES ------------------------------------------------------
my $file_dat = $dir_data."meso_seq_calls_matrix.txt";
my $file_prb = $dir_data."meso_pem_cnv_probe_medians_call.tsv";
my $file_seg = $dir_data."meso_pem_cnv_seg_mean_call_matrix.tsv";

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
	#my $ary_size = 53906;
	
	#print $ary_size, "\n";
	
	my $sample_name = $fields[0];
	
	print FILE_PRB $sample_name;
	print FILE_SEG $sample_name;
	
	### Probe median 
	for(my $i = 1; $i <= $ary_size; $i = $i + 2){
		print FILE_PRB "\t", $fields[$i];
	}
	print FILE_PRB "\n";
	
	### Segment values
	for(my $i = 2; $i <= $ary_size; $i = $i + 2){
		print FILE_SEG "\t", $fields[$i];
	}
	print FILE_SEG "\n";
}
close FILE_DAT;
close FILE_PRB;
close FILE_SEG;

exit;


#!/usr/bin/perl

use strict;
use warnings;

### SET GLOBAL PATHS -------------------------------------------------
my $dir_work = "/collinsgroup/Raunak/projects/MESO_peritoneal/data/cnv/seq_call_refseq_genes_meso/";

### DEFINE FILES -----------------------------------------------------
my $file_output = $dir_work."cnv_status_pem_nexus_export_transposed.tsv";

### READ FILE ---------------------------------------------------------
#open(FILE_DAT, $file_data) or die "Error: Cannot open data file";
open(FILE_OUT, ">$file_output") or die "Error: Cannot open output data file";


my(%data);          # main storage
my($maxcol) = 0;
my($rownum) = 0;
while (<>){
    #my(@row) = split /\s+/;
    my(@row) = split /\t/;
    my($colnum) = 0;
    foreach my $val (@row){
        $data{$rownum}{$colnum++} = $val;
    }
    $rownum++;
    $maxcol = $colnum if $colnum > $maxcol;
}

my $maxrow = $rownum;
for (my $col = 0; $col < $maxcol; $col++){
    for (my $row = 0; $row < $maxrow; $row++){
        printf FILE_OUT "%s%s", ($row == 0) ? "" : "\t",
        #printf "%s%s", ($row == 0) ? "" : "\t",
                defined $data{$row}{$col} ? $data{$row}{$col} : "";
    }
    print FILE_OUT "\n";
    #print "\n";
}

exit;


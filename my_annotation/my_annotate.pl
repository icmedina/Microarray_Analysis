#!/usr/bin/perl

# Program that annotate a list of selected genes based on their rpv
# using probe id and annotation metadata
# Isidro C. Medina Jr., 2008


open(RPV, "my_annotation/temp/out.csv"); @rpv = <RPV>; close(RPV);
open(ID, "my_annotation/metadata/probe_id.txt"); @probe_id = <ID>; close(ID);

# Get the probe id of the top 100 genes
open(ANNOTATE, ">my_annotation/temp/out.txt");
for ($n=0; $n<=99;$n++){
  $line = $rpv[$n];
 ($gene_index,$rpv) = split(/,/,$line);

    foreach $line_probe(@probe_id) {
      chomp($gene_index,$num,$probe,$rpv);
      ($num,$probe) = split(/\t/,$line_probe);

      if ($gene_index eq $num){
       print ANNOTATE "$gene_index\t$rpv\t$probe";
      }
    }
} close(ANNOTATE);

open(ANNOTATED, "my_annotation/temp/out.txt"); @annotated1 = <ANNOTATED>;
  close(ANNOTATED);
open(ANN, "my_annotation/metadata/annotation.txt"); @annotation = <ANN>; close(ANN);

open(ANNOTATION, ">my_annotation/biomarkers.xls");
 foreach $line1(@annotated1) {
 ($gene_index,$rpv,$probe) = split(/\t/,$line1);

    foreach $line2(@annotation) {
      chomp($gene_index,$rpv,$probe,$Affy_ProbeSet,$GeneName,$UnigeneComment,$Unigene,$LocusLink,$Chromosome,$CytoBand);
      ($Affy_ProbeSet,$GeneName,$UnigeneComment,$Unigene,$LocusLink,$Chromosome,$CytoBand) = split(/\t/,$line2);

      if ($probe eq $Affy_ProbeSet){
      print ANNOTATION "$rpv\t$gene_index\t$GeneName\t$UnigeneComment\t$LocusLink\n";
      }
    }
 } close(ANNOTATION);
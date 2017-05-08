#!/usr/bin/perl

open(RPV, "my_annotation/input/coef.txt"); @rpv = <RPV>; close(RPV);

# Get the probe id of the top 100 genes
open(ANNOTATE, ">my_annotation/ice.txt");
    foreach $line_probe(@rpv) {
      chomp($line_probe);
          
       print ANNOTATE "$line_probe, ";
      
    }
 close(ANNOTATE);


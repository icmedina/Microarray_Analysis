#!/usr/bin/perl

$input_rpv ="pRPV.csv";
$input_coef = "coef.txt";      # file that contains discriminant coefficients

$output = "coef-lambda.txt";
&get_index;


sub get_index{
  open(COEFF, "my_annotation/input/coef.txt"); @coef = <COEFF>; close(COEFF);
  open(RPV, "my_annotation/input/$input_rpv"); @rpv = <RPV>; close(RPV);

  open(DCOEF, ">my_annotation/$output");


    foreach $line_coef(@coef) {
    $line_coef =~ s/V//;

      chomp($line_coef);
      ($var,$coefficient) = split(/\t/,$line_coef);
       $index = $var - 2;

       foreach $line_rpv(@rpv){

 ($gene_index,$rpv) = split(/,/,$line_rpv);

      if ($index eq $gene_index){
       print DCOEF "$gene_index\t$index\t$coefficient\n";
      }
    }
} close(DCOEF);}
#!/usr/bin/perl

$input_rpv ="nRPV.csv";
$input_coef = "coef.txt";      # file that contains discriminant coefficients
$input_wilks = "wilks.txt";   # file that conatins wilks' lambda

$output = "coef-lambda.txt";

&get_wilks;
&get_index;

sub get_wilks{
  open(COEFF, "my_annotation/input/$input_coef"); @coef = <COEFF>; close(COEFF);
  open(WILKS, "my_annotation/input/$input_wilks"); @wilks = <WILKS>; close(WILKS);

  open(DCOEF, ">my_annotation/temp/var-wilks.txt");
foreach $line(@wilks){
     $line =~ s/ //g;
    chomp ($line);
   ($var_wilk,$tolerance,$F2remove,$wilks) = split(/\t/,$line);
    
    foreach $line_coef(@coef) {
    $line_coef =~ s/ //g;

      chomp($line_coef);
      ($var,$coefficient) = split(/\t/,$line_coef);
              
      if ($var eq $var_wilk){
      print DCOEF "$var\t$coefficient\t$tolerance\t$F2remove\t$wilks\n";
      }
    }
} close(DCOEF);}

sub get_index{
  open(COEFF, "my_annotation/temp/var-wilks.txt"); @coef = <COEFF>; close(COEFF);
  open(RPV, "my_annotation/input/$input_rpv"); @rpv = <RPV>; close(RPV);

  open(DCOEF, ">my_annotation/$output");
for ($n=0; $n<=100;$n++){
  $line = $rpv[$n];
 ($gene_index,$rpv) = split(/,/,$line);

    foreach $line_coef(@coef) {
    $line_coef =~ s/V//;
      
      chomp($line_coef);
      ($var,$coefficient,$tolerance,$F2remove,$wilks) = split(/\t/,$line_coef);
       $index = $var - 2;
       
      if ($index eq $n){
       print DCOEF "$gene_index\t$coefficient\t$tolerance\t$F2remove\t$wilks\n";
      }
    }
} close(DCOEF);}
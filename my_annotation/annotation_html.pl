#!/usr/bin/perl

# A program that output an HTML with annotation and links to NCBI
# Isidro C. Medina Jr., 2008
$condition = "0";

chomp($condition);
if ($condition eq "0"){
  $st = "n"; $state = "Down";}
elsif ($condition eq "1")
  {$st = "p"; $state = "Up";}

# Specify Input, output and Task Name
$input_coef = "coef.txt";
$input_rpv = "$st"."RPV.csv";
$task_name = "Serous Ovarian Cancer: $state-regulated Genes";
$out_path = "my_annotation/";
$annotated_file = "annotated-"."$state.txt";
$output_html = "annotation-"."$state.html";

$index_coef = "index-coef.txt";  # temporary file containing Gene Index and Coefficients
$probe_ids = "probe_id.txt"; 
$annotation_data = "annotation.txt";

&get_coef; 
&header; &body; &foot;
&get_annotation;

# GET GENE INDEX AND STANDARDIZED DISCRIMINANT COEFFICIENTS
sub get_coef {
open(COEFF, "my_annotation/input/$input_coef"); @coef = <COEFF>; close(COEFF);
open(RPV, "my_annotation/input/$input_rpv"); @rpv = <RPV>; close(RPV);

open(DCOEF, ">my_annotation/temp/$index_coef");
for ($n=0; $n<=100;$n++){
  $line = $rpv[$n];
 ($gene_index,$rpv) = split(/,/,$line);

    foreach $line_coef(@coef) {
    $line_coef =~ s/V//;
    $line_coef =~ s/ //g;
    
      ($num,$probe) = split(/\t/,$line_probe);
      chomp($line_coef);
   
      ($var,$coefficient) = split(/\t/,$line_coef);
       $var = $var - 2;
       
      if ($var eq $n){
       print DCOEF "$gene_index\t$coefficient\n";
      }
    }
} close(DCOEF);}

# GENERATE HTML REPORT
sub header{
open(ANNOTATE, ">$out_path$output_html");
  print ANNOTATE <<HEADER;
    <html><head></head><body bgcolor="#eeffff">
    <center>
    <h2>&nbsp;</h2>
    <h2>Biomarker Annotation Report</h2>
    <h2>$task_name</h2>
    <hr></center>
    <div align="center">
    <table border="1" cellpadding="3" cellspacing="0">
    <tbody>
      <tr> 
        <td colspan="3" align="center" bgcolor="#ffeeff">Feature selection method: 
          <b> Wavelet Variance</b></td>
      </tr>
      <tr align="center" bgcolor="#eeeeff"> 
        <td><strong>Gene Symbol<br>
          </strong></td>
        <td><strong>UniGene Description</strong></td>
        <td><strong>Locus Link</strong></td>
      </tr>
HEADER
close(ANNOTATE);}

sub body{
open(RPV, "my_annotation/temp/$index_coef"); @rpv = <RPV>; close(RPV);
open(ID, "my_annotation/metadata/$probe_ids"); @probe_id = <ID>; close(ID);

# Get the probe id of the top 100 genes
open(ANNOTATION, ">my_annotation/temp/out.txt");
foreach $line(@rpv){
 ($gene_index,$rpv) = split(/\t/,$line);

    foreach $line_probe(@probe_id) {
      chomp($gene_index,$num,$probe,$rpv);
      ($num,$probe) = split(/\t/,$line_probe);

      if ($gene_index eq $num){
       print ANNOTATION "$gene_index\t$rpv\t$probe";
      }
    }
} close(ANNOTATION);

open(ANNOTATED, "my_annotation/temp/out.txt"); @annotated1 = <ANNOTATED>;
  close(ANNOTATED);
open(ANN, "my_annotation/metadata/$annotation_data"); @annotation = <ANN>; close(ANN);

open(ANNOTATE, ">>$out_path$output_html");
 foreach $line1(@annotated1) {
 ($gene_index,$rpv,$probe) = split(/\t/,$line1);

    foreach $line2(@annotation) {
      chomp($gene_index,$rpv,$probe,$Affy_ProbeSet,$GeneName,$UnigeneComment,$Unigene,$LocusLink,$Chromosome,$CytoBand);
      ($Affy_ProbeSet,$GeneName,$UnigeneComment,$Unigene,$LocusLink,$Chromosome,$CytoBand) = split(/\t/,$line2);
      
      if (($GeneName ne "NULL") && ($probe eq $Affy_ProbeSet)){
        print ANNOTATE<<BODY;
          <tr align="center" bgcolor="#eeeeee">
            <td align="left"><div align="center">
           <a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=Nucleotide&amp;cmd=search&amp;term=$GeneName&amp;tool=gquery"> 
              $GeneName</a></div></td>
           
            <td align="left">$UnigeneComment</td>
            <td align="left"><div align="center">
            <a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=Gene&amp;cmd=search&amp;term=$LocusLink&amp;tool=gquery"> 
              $LocusLink</a></div></td>
          </tr>
BODY
      }
    }
 }
close(ANNOTATE);}

sub foot{
open(ANNOTATE, ">>$out_path$output_html");
  print ANNOTATE<<FOOTER;
    </tbody>
    </table>
    </div>
    <hr>
    <div align="center"><br>
    <font size="-1">COPYRIGHT, 2008</font><br>
    <strong><font size="3">Isidro C. Medina Jr</font></strong><font size="3">. 
    </font></div>
    </body></html>
FOOTER
close(ANNOTATE);}

# OUTPUT REPORT TO A TEXT FILE
sub get_annotation{
open(IDCOEF, "my_annotation/temp/$index_coef"); @id_coef = <IDCOEF>; close(IDCOEF);
open(ID, "my_annotation/metadata/$probe_ids"); @probe_id = <ID>; close(ID);

# Get the probe id of the top 100 genes
open(ANNOTATE, ">my_annotation/temp/out.txt");
foreach $line(@id_coef){
 ($gene_index,$coef) = split(/\t/,$line);

    foreach $line_probe(@probe_id) {
      chomp($gene_index,$num,$probe,$coef);
      ($num,$probe) = split(/\t/,$line_probe);

      if ($gene_index eq $num){
       print ANNOTATE "$gene_index\t$coef\t$probe";
      }
    }
} close(ANNOTATE);

open(ANNOTATED, "my_annotation/temp/out.txt"); @annotated1 = <ANNOTATED>;
  close(ANNOTATED);
open(ANN, "my_annotation/metadata/$annotation_data"); @annotation = <ANN>; close(ANN);

open(ANNOTATION, ">my_annotation/$annotated_file");
 foreach $line1(@annotated1) {
 ($gene_index,$id_coef,$probe) = split(/\t/,$line1);

    foreach $line2(@annotation) {
      chomp($gene_index,$id_coef,$probe,$Affy_ProbeSet,$GeneName,$UnigeneComment,$Unigene,$LocusLink,$Chromosome,$CytoBand);
      ($Affy_ProbeSet,$GeneName,$UnigeneComment,$Unigene,$LocusLink,$Chromosome,$CytoBand) = split(/\t/,$line2);

      if ($probe eq $Affy_ProbeSet){
      print ANNOTATION "$gene_index\t$id_coef\t$GeneName\t$UnigeneComment\t$LocusLink\n";
      }
    }
 } close(ANNOTATION);}
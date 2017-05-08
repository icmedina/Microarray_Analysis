#!D:\Perl64\bin\perl.exe
# Program that removes suffixes in the accession number of Affymetrix ID
# Isidro C. Medina Jr., 2008

open (TXTFILE, "probe_id_in.txt");   #input file
@lines = <TXTFILE>;
close(TXTFILE);

open (OUTPUT, ">>probe_id_out.txt");        #output file
  foreach $line (@lines) {
    $line=~s/\_.*//;
    print OUTPUT "$line";}
close(OUTPUT);
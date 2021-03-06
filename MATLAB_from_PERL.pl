#!D:\perl\bin\perl.exe -w
use Win32::OLE;
use Win32::OLE::Variant;

# Perl script to execute commands in Matlab.

# Use existing instance if Matlab is already running.
eval {$ml = Win32::OLE->GetActiveObject('Matlab.Application')};
die "Matlab not installed" if $@;
unless (defined $ml) {
   $ml = Win32::OLE->new('Matlab.Application')
      or die "Oops, cannot start Matlab";
}

# Executing MATLAB commands
# MATLAB and return the status.
@exe_commands = ("rma = dlmread(\'RMA.csv\', \',\');",
             "rma");

# Send the commands to Matlab
foreach $exe_command (@exe_commands)
{  $status = &send_to_matlab('Execute', $exe_command);
   print "Matlab status = ", $status, "\n";
}

sub send_to_matlab
{  my ($call, @command) = @_;
   my $status = 0;
   print "\n>> $call( @command )\n";
   $result = $ml->Invoke($call, @command);
   if (defined($result))
   {   unless ($result =~ s/^.\?{3}/Error:/)
       {  print "$result\n" unless ($result eq "");
       }
       else
       {  print "$result\n";
          $status = -1;
       }
   }
   return $status;
}

# Examples of passing variables between MATLAB and Perl.
#
# MATLAB supoprts passing character arrays directly with the following syntax:
#
# PutCharArray([in] BSTR name, [in] BSTR workspace, [in] BSTR string);
# GetCharArray([in] BSTR name, [in] BSTR workspace, [out] BSTR string);

&send_to_matlab('PutCharArray', 'centralDogma', 'base', 'DNA->RNA->Protein.');
&send_to_matlab('GetCharArray', 'centralDogma', 'base');

# Numeric arrays can be passed by reference in a SAFEARRAY using the
# PutFullMatrix and GetFullMatrix functions.
#
# PutFullMatrix([in] BSTR name, [in] BSTR workspace, [in] BSTR data);
# GetFullMatrix([in] BSTR varname, [in] BSTR workspace, [out] BSTR retdata);

$mReal = Variant(VT_ARRAY|VT_R8, 4, 4);
$mImag = Variant(VT_ARRAY|VT_R8, 4, 4);

$mReal->Put([[0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]]);
print "\n>> PutFullMatrix( 'magicArray', 'base', ",'$mReal, $mImag'," )\n";
$ml->PutFullMatrix('magicArray', 'base', $mReal, $mImag);
$ml->Execute('magicArray = magic(4)');

$m2Real = Variant(VT_ARRAY|VT_R8|VT_BYREF,4,4);
$m2Imag = Variant(VT_ARRAY|VT_R8|VT_BYREF,4,4);
print "\n>> GetFullMatrix( 'magicArray', 'base', ",'$m2Real, $m2Imag'," )\n";
$ml->GetFullMatrix('magicArray', 'base', $m2Real, $m2Imag);

for ($i = 0; $i < 4; $i++) {
 printf "%3d %3d %3d %3d\n", $m2Real->Get($i,0), $m2Real->Get($i,1),
                                $m2Real->Get($i,2), $m2Real->Get($i,3);
}

# Additionally, you can use Variants to send scalar variables by reference
# to MATLAB for all data types except sparse arrays and function handles through
# PutWorkspaceData:
# PutWorkspaceData([in] BSTR name, [in] BSTR workspace, [in] BSTR data);
#
# Results are passed back to Perl directly with GetVariable:
# HRESULT = GetVariable([in] BSTR Name, [in] BSTR Workspace);

# Create and initialize a date Variant.
$dnaDate = Variant->new(VT_DATE|VT_BYREF, 'Feb 28, 1953');

&send_to_matlab('PutWorkspaceData', 'dnaDate', 'base', $dnaDate);
&send_to_matlab('Execute', 'dnaDate');

# Create and initialize a new string Variant.
$aminoString = Variant->new(VT_BSTR|VT_BYREF, 'matlap');
&send_to_matlab('PutWorkspaceData', 'aminoAcids', 'base', $aminoString);

# Change the value in MATLAB
&send_to_matlab('Execute', "aminoAcids = 'ARNDCQEGHILKMFPSTWYV';");

# Bring the new value back
$aa = $ml->GetVariable('aminoAcids', 'base');
printf "Amino acid codes: %s\n", $aa;

undef $ml; # close Matlab if we opened it
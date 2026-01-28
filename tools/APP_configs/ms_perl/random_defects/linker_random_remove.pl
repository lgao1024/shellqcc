#!/usr/bin/perl
use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

my $doc = $Documents{"UiO-66-SP_5.xsd"};
my $NumberOfDeleted=0; #input: number of deleted ligands.
my $bonds = $doc->UnitCell->Bonds;
our @fragments=();
our $fragmentsnum=0;

foreach my $bond (@$bonds) {
    if (($bond->Atom1->ElementSymbol eq "N") and ($bond->Atom2->ElementSymbol eq "N")) {
        my $fragmentAtoms=$bond->Atom1->FragmentUnitCell->Atoms;
        $fragmentsnum=$fragmentsnum+1;
        push(@fragments, [$fragmentAtoms,$fragmentsnum]);
    }
}

my $range = $fragmentsnum;
# building deleting list
our @randnumlist=();
our $randnumconter=0;
do{
    my $random_number = int(rand($range));
    my $flag=0;
    foreach my $j (@randnumlist){
        if($random_number eq $j){$flag=1;}
    }
    if($flag eq "0"){
        push(@randnumlist,$random_number);
        $randnumconter=$randnumconter+1;
    }
}while( $randnumconter< $NumberOfDeleted );

#check the length of randnumlist
my $ii=0;
foreach my $i (@randnumlist) {
    $ii=$ii+1;
}
if ($ii eq $NumberOfDeleted){
    #delete the ligand
    foreach my $i (@randnumlist) {
        @fragments[$i]->[0]->Delete;
        print "Deleting ligand NO.:",$i,"\n";
    }
}
print "fragments number is:###",$fragmentsnum,"###","\n";
print "% of ligands had been deleted:",100*$NumberOfDeleted/$fragmentsnum,"%","\n";
#calculate bonds
#$doc->CalculateBonds;
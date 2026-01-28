#!/usr/bin/perl
use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

my $doc = $Documents{"UiO-66-SP_5.xsd"};
my $NumberOfDeleted = 0;  # Input: number of Zr atoms to delete

# Step 1: Collect all Zr atoms
my $atoms = $doc->UnitCell->Atoms;
our @zr_atoms = ();
our $zr_count = 0;

foreach my $atom (@$atoms) {
    if ($atom->ElementSymbol eq "Zr") {
        push(@zr_atoms, $atom);
        $zr_count++;
    }
}

# Step 2: Randomly select Zr atoms to delete
our @rand_indices = ();
our $counter = 0;

do {
    my $random_index = int(rand($zr_count));
    my $already_selected = 0;
    foreach my $j (@rand_indices) {
        if ($random_index eq $j) {
            $already_selected = 1;
            last;
        }
    }
    if ($already_selected == 0) {
        push(@rand_indices, $random_index);
        $counter++;
    }
} while ($counter < $NumberOfDeleted);

# Step 3: Delete selected Zr atoms
foreach my $i (@rand_indices) {
    $zr_atoms[$i]->Delete;
    print "Deleted Zr atom index: $i\n";
}

print "Total Zr atoms before deletion: $zr_count\n";
print "Deleted $NumberOfDeleted Zr atoms (", 100 * $NumberOfDeleted / $zr_count, "%)\n";

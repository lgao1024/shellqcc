#!/usr/bin/perl
use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

my $doc = $Documents{"SC6_linker_468-accessibleCu-label.xsd"};
my $NumberOfReplaced = 33;  # Input: number of Zr atoms to replace with Cu

# Step 0: basic sanity check on the document
$doc or die "Failed to open document.\n";

# Step 1: Collect all Zr atoms
my $atoms = $doc->UnitCell->Atoms;  # array ref of Atom objects
my @zr_atoms;
my $zr_count = 0;

foreach my $atom (@$atoms) {
    if ($atom->ElementSymbol eq 'Zr') {
        push @zr_atoms, $atom;
        $zr_count++;
    }
}

die "Requested $NumberOfReplaced replacements, but only $zr_count Zr atoms found\n"
    if $NumberOfReplaced > $zr_count;

# Step 2: Randomly select Zr atoms to replace (no duplicates)
my @rand_indices;
my $counter = 0;

while ($counter < $NumberOfReplaced) {
    my $random_index = int(rand($zr_count));
    my $already_selected = 0;
    foreach my $j (@rand_indices) {
        if ($random_index == $j) {  # numeric comparison
            $already_selected = 1;
            last;
        }
    }
    if (!$already_selected) {
        push @rand_indices, $random_index;
        $counter++;
    }
}

# Step 3: Replace selected Zr atoms with Cu
foreach my $i (@rand_indices) {
    # In MaterialsScript, ElementSymbol is a property (getter with 0 args). 
    # Use assignment to set it, not a method call with a parameter.
    $zr_atoms[$i]->ElementSymbol = 'Cu';
    print "Replaced Zr atom index $i with Cu\n";
}

print "Total Zr atoms before replacement: $zr_count\n";
printf "Replaced %d Zr atoms with Cu (%.2f%%)\n",
        $NumberOfReplaced, 100 * $NumberOfReplaced / $zr_count;

# Optional: save to a new file
my $out = "SC6_linker_468-accessibleCu-label_REPLACED.xsd";
$doc->SaveAs($out);
print "Saved: $out\n";

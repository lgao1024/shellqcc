#!perl
use strict;
use warnings;
use MaterialsScript qw(:all);

# =======================
# Input documents
# =======================
my $struct = $Documents{"UiO-66.xsd"};     # 3D Atomistic Document
my $exp    = $Documents{"UiO-66pw.xcd"};   # Chart Document (powder pattern)

# =======================
# Helpers
# =======================
sub get_prop {
    my ($obj, @names) = @_;
    for my $n (@names) {
        my $v = eval { $obj->$n };
        return $v if defined $v;
    }
    return undef;
}

sub fmt6 {
    my ($v) = @_;
    return "(n/a)" if !defined $v;
    return sprintf("%.6f", $v);
}

# Convert a KeyValuePairs-like object to a Perl hash (works across MS builds)
sub kvp_to_hash {
    my ($kvp) = @_;
    my %h;
    return %h unless $kvp;

    # Pattern A: Keys() -> arrayref
    my $keys = eval { $kvp->Keys };
    if ($keys && ref($keys) eq "ARRAY") {
        foreach my $k (@$keys) {
            my $v = eval { $kvp->Value($k) };
            $v = eval { $kvp->Item($k) } unless defined $v;
            $v = eval { $kvp->Get($k) }  unless defined $v;
            $h{$k} = $v if defined $v;
        }
        return %h if keys %h;
    }

    # Pattern B: Count + Key(i) + Value(i)
    my $count = eval { $kvp->Count };
    if (defined $count && $count =~ /^\d+$/) {
        for (my $i = 0; $i < $count; $i++) {
            my $k = eval { $kvp->Key($i) };
            $k = eval { $kvp->Keys($i) } unless defined $k;
            my $v = eval { $kvp->Value($i) };
            $v = eval { $kvp->Item($i) } unless defined $v;
            $h{$k} = $v if defined $k && defined $v;
        }
        return %h if keys %h;
    }

    # Pattern C: AsHash / Hash / ToHash
    my $href = eval { $kvp->AsHash };
    $href = eval { $kvp->Hash }   unless ($href && ref($href) eq "HASH");
    $href = eval { $kvp->ToHash } unless ($href && ref($href) eq "HASH");
    if ($href && ref($href) eq "HASH") {
        %h = %$href;
        return %h if keys %h;
    }

    return %h;
}

# =======================
# Reflex: enable lattice refinement
# =======================
my $prep = Modules->Reflex->Preparation;
$prep->SetRefineLattice($struct, "Yes");

# =======================
# Reflex settings
# =======================
Modules->Reflex->ChangeSettings([
  WriteLevel                => "Verbose",

  RefinementMethod          => "Pawley",
  InstrumentGeometry        => "Bragg-Brentano",
  ProfileFunction           => "Pseudo-Voigt",

  BackgroundOrder           => 20,
  ShowBackground            => "Yes",
  RefineBackground          => "Yes",

  TwoThetaMin               => 4,
  TwoThetaMax               => 40,

  FixFractionalCoordinates  => "Yes",

  NumRefinementCycles       => 10,
  NumFunctionEvaluations    => 30,

  #ZeroPoint                 => 0.0,
  #ProfileU                  => 0.01,
  #ProfileV                  => -0.001,
  #ProfileW                  => 0.002,
  #ProfileNA                 => 0.5,
  
  ZeroPoint                 => "-0.117208",
  ProfileU                  => "0.410190",
  ProfileV                  => "-0.219287",
  ProfileW                  => "0.037100",
  ProfileNA                 => "0.681390",

  RefineZeroPoint           => "Yes",
  RefineProfileU            => "Yes",
  RefineProfileV            => "Yes",
  RefineProfileW            => "Yes",
  RefineProfileNA           => "Yes",

  ShowDifferencePattern     => "Yes",
]);

# =======================
# Run refinement
# =======================
my $results = Modules->Reflex->PowderRefinement->Run($struct, $exp);

# =======================
# Output Rp / Rwp
# =======================
printf "Rp  = %.4f %%\n",  $results->Rp;
printf "Rwp = %.4f %%\n\n", $results->Rwp;

# =======================
# Output refined lattice (numeric)
# =======================
my $refStruct = eval { $results->Structure } || $struct;
my $lat3d     = eval { $refStruct->Lattice3D };

print "Refined lattice parameters (numeric):\n";
if ($lat3d) {
    my $a  = get_prop($lat3d, qw(LengthA A));
    my $b  = get_prop($lat3d, qw(LengthB B));
    my $c  = get_prop($lat3d, qw(LengthC C));
    my $al = get_prop($lat3d, qw(AngleAlpha Alpha));
    my $be = get_prop($lat3d, qw(AngleBeta Beta));
    my $ga = get_prop($lat3d, qw(AngleGamma Gamma));

    print "  a     = " . fmt6($a)  . "\n";
    print "  b     = " . fmt6($b)  . "\n";
    print "  c     = " . fmt6($c)  . "\n";
    print "  alpha = " . fmt6($al) . "\n";
    print "  beta  = " . fmt6($be) . "\n";
    print "  gamma = " . fmt6($ga) . "\n";
} else {
    print "  (no Lattice3D found on structure)\n";
}

# =======================
# Output refined profile parameters (copy/paste, using RefinedValues)
# =======================
print "\nRefined profile parameters (copy/paste):\n";

my $rv = eval { $results->RefinedValues };
my %rvh = kvp_to_hash($rv);

# In your build, these keys exist exactly as follows:
# ZeroPoint, ProfileU, ProfileV, ProfileW, ProfileNA
my $zp = $rvh{"ZeroPoint"};
my $u  = $rvh{"ProfileU"};
my $v  = $rvh{"ProfileV"};
my $w  = $rvh{"ProfileW"};
my $na = $rvh{"ProfileNA"};

# Print exactly in the format you requested: => "0.01",
printf "  ZeroPoint                 => \"%.6f\",\n", $zp if defined $zp;
printf "  ProfileU                  => \"%.6f\",\n", $u  if defined $u;
printf "  ProfileV                  => \"%.6f\",\n", $v  if defined $v;
printf "  ProfileW                  => \"%.6f\",\n", $w  if defined $w;
printf "  ProfileNA                 => \"%.6f\",\n", $na if defined $na;

if (!defined $zp && !defined $u && !defined $v && !defined $w && !defined $na) {
    print "  (not found in RefinedValues; check that RefineZeroPoint/RefineProfile* are Yes)\n";
}

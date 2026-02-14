#!perl

use strict;
use Scalar::Util qw(refaddr);
use POSIX qw(floor);
use MaterialsScript qw(:all);

$| = 1;

# ============================================================
# USER-TUNABLE PARAMETERS
# ============================================================

# Manual fallback input. Leave empty to rely on ActiveDocument.
my $INPUT = "";

# Prefer active document when script is run via User menu command ("Active document only")
# 1 = ON (default): try Documents->ActiveDocument first, fallback to $INPUT if not empty
# 0 = OFF: always use $INPUT
my $USE_ACTIVE_DOCUMENT = 1;

# --- Disconnect options BEFORE counting fragments ---
my $CUT_METAL_O = 1;
my $CUT_METAL_N = 1;
my $N_DEG_MAX_FOR_CUT = 2;

# --- Ligand-side O identification (only for CUT_METAL_O=1) ---
my $MARKER = "Xe";
my $USE_GEOM_FALLBACK = 1;
my $GEOM_CUTOFF_A = 1.90;

# --- Defect (random deletion) settings ---
my $DEFECT_LIGAND_FRAC  = 0.0;
my $DEFECT_CLUSTER_FRAC = 0.0;

# Random seed:
# -1 => time-based, integer => reproducible
my $RANDOM_SEED = 12345;

# --- Post-processing ---
my $RECALC_BONDS_AT_END = 0;

# --- Output options (ONLY affects printing; does NOT change logic) ---
my $REPORT_WIDTH = 72;
my $VERBOSE = 0;
my $SHOW_FRAGMENT_ID_LIST = 0;
my $MAX_ID_LIST = 80;
my $REPORT_STYLE = "pretty";  # "pretty" or "tsv"

# --- Metal element table (broad default) ---
my @METAL_TABLE = qw(
    Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In
    Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi
    U Th Pa Np Pu Am Cm Bk Cf Es Fm Md No Lr
);

# ============================================================
# END USER PARAMETERS
# ============================================================

# -------------------------
# Helper: metal lookup
# -------------------------
my %METAL = map { $_ => 1 } @METAL_TABLE;
sub is_metal_el { my ($el) = @_; return $METAL{$el} ? 1 : 0; }

# -------------------------
# Output helpers (ONLY printing)
# -------------------------
sub hr {
    my ($ch) = @_;
    $ch = "=" unless defined $ch && length($ch) == 1;
    print $ch x $REPORT_WIDTH, "\n";
}

sub section {
    my ($title) = @_;
    if ($REPORT_STYLE eq "tsv") {
        print "# $title\n";
        return;
    }
    print "\n";
    hr("=");
    print "$title\n";
    hr("=");
}

sub kv {
    my ($k, $v) = @_;
    $v = "" unless defined $v;
    if ($REPORT_STYLE eq "tsv") {
        $k =~ s/\t/ /g;
        $v =~ s/\t/ /g;
        print "$k\t$v\n";
        return;
    }
    printf("  %-28s : %s\n", $k, $v);
}

sub pct_str {
    my ($n, $d) = @_;
    $d = 0 + $d;
    return "0.00%" if $d <= 0;
    my $p = 100.0 * $n / $d;
    return sprintf("%.2f%%", $p);
}

sub list_preview {
    my ($label, $arr_ref) = @_;
    return unless $SHOW_FRAGMENT_ID_LIST;
    my @ids = @$arr_ref;
    my $n = scalar(@ids);
    my $show = $n;
    $show = $MAX_ID_LIST if $show > $MAX_ID_LIST;

    my @head = ();
    @head = @ids[0 .. ($show-1)] if $show > 0;
    my $tail = ($n > $show) ? " ... (+".($n-$show).")" : "";

    kv($label, ($show ? join(",", @head).$tail : "(none)"));
}

# -------------------------
# Helpers: stable atom key in periodic systems
# -------------------------
sub mod1 {
    my ($x) = @_;
    $x = $x - floor($x);
    $x = 0.0 if abs($x - 1.0) < 1e-9;
    return $x;
}

sub atom_key {
    my ($a) = @_;
    my $e = $a->ElementSymbol;
    my $p = $a->FractionalXYZ;
    return sprintf("%s|%.6f|%.6f|%.6f", $e, mod1($p->X), mod1($p->Y), mod1($p->Z));
}

# -------------------------
# Helper: shuffle array in-place (Fisher–Yates)
# -------------------------
sub shuffle_in_place {
    my ($arr_ref) = @_;
    for (my $i = @$arr_ref - 1; $i > 0; $i--) {
        my $j = int(rand($i + 1));
        ($arr_ref->[$i], $arr_ref->[$j]) = ($arr_ref->[$j], $arr_ref->[$i]);
    }
}

sub clamp01 {
    my ($x) = @_;
    $x = 0 if $x < 0;
    $x = 1 if $x > 1;
    return $x;
}

# ============================================================
# 0) Resolve document source (ActiveDocument preferred)
#    (match MS example style for User Script activation)
# ============================================================
section("MOF Fragment Report");
kv("INPUT (manual fallback)", ($INPUT eq "" ? "(empty)" : $INPUT));
kv("USE_ACTIVE_DOCUMENT", $USE_ACTIVE_DOCUMENT ? "ON" : "OFF");

my ($doc, $doc_source) = (undef, "");

if ($USE_ACTIVE_DOCUMENT) {
    eval {
        $doc = Documents->ActiveDocument;
    };
    if ($@ || !$doc) {
        $doc = undef;
    } else {
        $doc_source = "ActiveDocument";
    }
}

if (!$doc) {
    if ($INPUT eq "") {
        die "Cannot open document: ActiveDocument unavailable and INPUT is empty.\n"
          . "Tip: run via User Script with 'Active document only', or set \$INPUT to a filename.\n";
    }
    $doc = $Documents{$INPUT};
    $doc_source = "INPUT";
}

die "Cannot open document (source=$doc_source).\n" unless $doc;

kv("Document source", $doc_source);
eval { kv("Document name", $doc->Name); };

# ============================================================
# 1) Open UnitCell
# ============================================================
my $uc = $doc->UnitCell;
die "No UnitCell found.\n" unless $uc;

my @atoms = @{ $uc->Atoms };
my @bonds = @{ $uc->Bonds };
die "No atoms found in UnitCell.\n" if @atoms == 0;

kv("CellFormula", $doc->SymmetrySystem->CellFormula);
kv("UnitCell atoms", scalar(@atoms));
kv("UnitCell bonds (initial)", scalar(@bonds));

# ============================================================
# A) DISCONNECT coordination links (optional)
# ============================================================

# Present elements
my %present;
foreach my $a (@atoms) { $present{ $a->ElementSymbol } = 1; }
my @present_els = sort keys %present;

# Metals present
my %present_metals;
foreach my $el (@present_els) { $present_metals{$el} = 1 if is_metal_el($el); }
my @metal_els = sort keys %present_metals;

section("Elements");
kv("Present elements", join(",", @present_els));
kv("Detected metal elements", (@metal_els ? join(",", @metal_els) : "(none)"));

# Ligand heavy elements (non-metal, non-H, non-O)
my %ligand_heavy;
foreach my $el (@present_els) {
    next if $el eq "H";
    next if $el eq "O";
    next if is_metal_el($el);
    $ligand_heavy{$el} = 1;
}
my @ligand_heavy_els = sort keys %ligand_heavy;

section("Disconnect Settings");
kv("CUT_METAL_O", ($CUT_METAL_O ? "ON" : "OFF"));
kv("CUT_METAL_N", ($CUT_METAL_N ? "ON" : "OFF"));
kv("MARKER", $MARKER);
kv("USE_GEOM_FALLBACK", ($USE_GEOM_FALLBACK ? "ON" : "OFF"));
kv("GEOM_CUTOFF_A (Å)", $GEOM_CUTOFF_A);
kv("N_DEG_MAX_FOR_CUT", $N_DEG_MAX_FOR_CUT);
kv("Ligand-heavy elements (for O marking)", (@ligand_heavy_els ? join(",", @ligand_heavy_els) : "(none)"));

# ---- A1) Mark ligand-side O as MARKER (only if CUT_METAL_O) ----
my %marked_O_by_addr;
my ($marked_from_bonds, $marked_from_geom, $ox_bonds_found) = (0, 0, 0);

if ($CUT_METAL_O) {
    die "Marker element '$MARKER' already exists; choose a different MARKER.\n" if $present{$MARKER};
    die "No ligand-heavy elements detected; cannot mark ligand-side O.\n"
        if (!@ligand_heavy_els && $USE_GEOM_FALLBACK);

    foreach my $bond (@bonds) {
        my $a1 = $bond->Atom1;
        my $a2 = $bond->Atom2;
        next unless $a1 && $a2;

        my $e1 = $a1->ElementSymbol;
        my $e2 = $a2->ElementSymbol;

        my $o;
        if ($e1 eq "O" && $ligand_heavy{$e2}) { $o = $a1; }
        elsif ($e2 eq "O" && $ligand_heavy{$e1}) { $o = $a2; }
        else { next; }

        $ox_bonds_found++;
        my $id = refaddr($o);
        next if $marked_O_by_addr{$id};

        $o->ElementSymbol = $MARKER;
        $marked_O_by_addr{$id} = 1;
        $marked_from_bonds++;
    }

    if ($USE_GEOM_FALLBACK) {
        my @X_atoms = grep { $ligand_heavy{ $_->ElementSymbol } } @atoms;
        my @O_atoms = grep { $_->ElementSymbol eq "O" } @atoms;

        my $r2 = $GEOM_CUTOFF_A * $GEOM_CUTOFF_A;
        my @Xx = map { $_->X } @X_atoms;
        my @Xy = map { $_->Y } @X_atoms;
        my @Xz = map { $_->Z } @X_atoms;

        for (my $i = 0; $i < @O_atoms; $i++) {
            my $o = $O_atoms[$i];
            my $id = refaddr($o);
            next if $marked_O_by_addr{$id};

            my ($ox, $oy, $oz) = ($o->X, $o->Y, $o->Z);
            my $hit = 0;

            for (my $j = 0; $j < @X_atoms; $j++) {
                my $dx = $ox - $Xx[$j];
                my $dy = $oy - $Xy[$j];
                my $dz = $oz - $Xz[$j];
                my $d2 = $dx*$dx + $dy*$dy + $dz*$dz;
                if ($d2 <= $r2) { $hit = 1; last; }
            }

            next unless $hit;
            $o->ElementSymbol = $MARKER;
            $marked_O_by_addr{$id} = 1;
            $marked_from_geom++;
        }
    }

    kv("O-(ligand heavy) bonds found", $ox_bonds_found);
    kv("Unique O marked (bond-based)", $marked_from_bonds);
    kv("Additional O marked (geom)", $marked_from_geom) if $USE_GEOM_FALLBACK;
    kv("Total unique marked O atoms", scalar(keys %marked_O_by_addr));
} else {
    kv("O marking", "SKIPPED (CUT_METAL_O=OFF)");
}

@bonds = @{ $uc->Bonds };

# ---- A2) Select bonds to delete: Metal–MARKER and/or Metal–N ----
my @to_delete;
foreach my $bond (@bonds) {
    my $a1 = $bond->Atom1;
    my $a2 = $bond->Atom2;
    next unless $a1 && $a2;

    my $e1 = $a1->ElementSymbol;
    my $e2 = $a2->ElementSymbol;

    if ($CUT_METAL_O) {
        if ( ($present_metals{$e1} && $e2 eq $MARKER) || ($present_metals{$e2} && $e1 eq $MARKER) ) {
            push @to_delete, $bond;
            next;
        }
    }

    if ($CUT_METAL_N) {
        my $n;
        if ($present_metals{$e1} && $e2 eq "N") { $n = $a2; }
        elsif ($present_metals{$e2} && $e1 eq "N") { $n = $a1; }
        else { next; }

        if ($N_DEG_MAX_FOR_CUT && $N_DEG_MAX_FOR_CUT > 0) {
            my $ndeg = scalar(@{ $n->Bonds });
            next if $ndeg > $N_DEG_MAX_FOR_CUT;
        }

        push @to_delete, $bond;
        next;
    }
}

kv("Bonds selected for deletion", scalar(@to_delete));

# ---- A3) Delete selected bonds ----
my ($deleted, $failed) = (0, 0);
foreach my $bond (@to_delete) {
    my $ok = 0;
    eval { $bond->Delete; $ok = 1; };
    if (!$ok) {
        eval { $uc->Bonds->Delete($bond); $ok = 1; };
    }
    $ok ? $deleted++ : $failed++;
}
kv("Bonds deleted", $deleted);
kv("Bonds failed to delete", $failed);

# ---- A4) Revert MARKER back to O ----
if ($CUT_METAL_O) {
    my $reverted = 0;
    foreach my $a (@{ $uc->Atoms }) {
        next unless $a->ElementSymbol eq $MARKER;
        $a->ElementSymbol = "O";
        $reverted++;
    }
    kv("Reverted marker atoms back to O", $reverted);
}

kv("UnitCell bonds (after disconnect)", scalar(@{ $uc->Bonds }));

# ============================================================
# B) BUILD fragment dictionaries by connected components (graph on UnitCell)
# ============================================================
section("Fragments");

my $atoms_uc = $uc->Atoms;

my %key2atom;
foreach my $a (@$atoms_uc) { $key2atom{ atom_key($a) } = $a; }

my %IS_METAL = map { $_ => 1 } @metal_els;

my %CLUSTER_FRAG = ();
my %LIGAND_FRAG  = ();
my @CLUSTER_IDS = ();
my @LIGAND_IDS  = ();

my $frag_total = 0;

my %visited;
my $next_cluster_id = 0;
my $next_ligand_id  = 0;

foreach my $a0 (@$atoms_uc) {
    my $k0 = atom_key($a0);
    next if $visited{$k0};

    my @queue = ($k0);
    $visited{$k0} = 1;

    my @frag_keys;
    while (@queue) {
        my $k = shift @queue;
        push @frag_keys, $k;

        my $a = $key2atom{$k};
        next unless defined $a;

        my $bonds_here = $a->Bonds;
        foreach my $bond (@$bonds_here) {
            my $x1 = $bond->Atom1;
            my $x2 = $bond->Atom2;

            my $k1 = atom_key($x1);
            my $k2 = atom_key($x2);
            my $kn = ($k1 eq $k) ? $k2 : $k1;

            next if $visited{$kn};
            next unless exists $key2atom{$kn};

            $visited{$kn} = 1;
            push @queue, $kn;
        }
    }

    # classify fragment (original logic preserved)
    my ($has_metal, $has_c) = (0, 0);
    foreach my $k (@frag_keys) {
        my $aa = $key2atom{$k};
        my $e  = $aa->ElementSymbol;
        $has_c = 1 if $e eq "C";
        $has_metal = 1 if $IS_METAL{$e};
        last if $has_metal;
    }

    $frag_total++;

    if ($has_metal) {
        my $id = $next_cluster_id++;
        $CLUSTER_FRAG{$id} = [ @frag_keys ];
        push @CLUSTER_IDS, $id;
    } elsif ($has_c) {
        my $id = $next_ligand_id++;
        $LIGAND_FRAG{$id} = [ @frag_keys ];
        push @LIGAND_IDS, $id;
    }
}

kv("Fragments total", $frag_total);
kv("Cluster fragments", scalar(@CLUSTER_IDS));
kv("Ligand fragments", scalar(@LIGAND_IDS));

# ============================================================
# C) RANDOM deletion of fragments by ratio (0..1)
# ============================================================
section("Defect Settings");

$DEFECT_LIGAND_FRAC  = clamp01($DEFECT_LIGAND_FRAC);
$DEFECT_CLUSTER_FRAC = clamp01($DEFECT_CLUSTER_FRAC);

if ($RANDOM_SEED < 0) {
    srand(time());
    kv("Random seed", "time-based");
} else {
    srand($RANDOM_SEED);
    kv("Random seed", $RANDOM_SEED);
}

my $n_lig_total = scalar(@LIGAND_IDS);
my $n_clu_total = scalar(@CLUSTER_IDS);

my $n_lig_del = int($DEFECT_LIGAND_FRAC  * $n_lig_total + 0.5);
my $n_clu_del = int($DEFECT_CLUSTER_FRAC * $n_clu_total + 0.5);

kv("Ligand delete plan",  sprintf("frac=%s -> %d / %d", $DEFECT_LIGAND_FRAC,  $n_lig_del, $n_lig_total));
kv("Cluster delete plan", sprintf("frac=%s -> %d / %d", $DEFECT_CLUSTER_FRAC, $n_clu_del, $n_clu_total));

my @lig_ids = @LIGAND_IDS;
my @clu_ids = @CLUSTER_IDS;

shuffle_in_place(\@lig_ids);
shuffle_in_place(\@clu_ids);

my @lig_del_ids = ($n_lig_del > 0) ? @lig_ids[ 0 .. ($n_lig_del-1) ] : ();
my @clu_del_ids = ($n_clu_del > 0) ? @clu_ids[ 0 .. ($n_clu_del-1) ] : ();

list_preview("Ligand fragment IDs to delete", \@lig_del_ids);
list_preview("Cluster fragment IDs to delete", \@clu_del_ids);

section("Deletion");

my $atoms_deleted = 0;
my $lig_frag_deleted_atoms = 0;
my $clu_frag_deleted_atoms = 0;

sub delete_fragment_atoms_by_keys {
    my ($frag_type, $frag_id, $keys_ref, $key2atom_ref) = @_;
    my $deleted_here = 0;

    foreach my $k (@$keys_ref) {
        my $a = $key2atom_ref->{$k};
        next unless defined $a;

        my $ok = 0;
        eval { $a->Delete; $ok = 1; };
        $deleted_here++ if $ok;
    }

    print "  [DELETE] $frag_type fragment id=$frag_id : atoms_deleted=$deleted_here\n" if $VERBOSE;
    return $deleted_here;
}

foreach my $id (@lig_del_ids) {
    my $d = delete_fragment_atoms_by_keys("LIGAND", $id, $LIGAND_FRAG{$id}, \%key2atom);
    $atoms_deleted += $d;
    $lig_frag_deleted_atoms += $d;
}

foreach my $id (@clu_del_ids) {
    my $d = delete_fragment_atoms_by_keys("CLUSTER", $id, $CLUSTER_FRAG{$id}, \%key2atom);
    $atoms_deleted += $d;
    $clu_frag_deleted_atoms += $d;
}

kv("Ligand frags deleted",  sprintf("%d / %d (%s)", scalar(@lig_del_ids), $n_lig_total, pct_str($n_lig_del, $n_lig_total)));
kv("Cluster frags deleted", sprintf("%d / %d (%s)", scalar(@clu_del_ids), $n_clu_total, pct_str($n_clu_del, $n_clu_total)));
kv("Atoms deleted (ligand)",  $lig_frag_deleted_atoms);
kv("Atoms deleted (cluster)", $clu_frag_deleted_atoms);
kv("Total atoms deleted", $atoms_deleted);

section("POST: RE-CALCULATE BONDS");
kv("RECALC_BONDS_AT_END", $RECALC_BONDS_AT_END ? "ON" : "OFF");

if ($RECALC_BONDS_AT_END) {
    my $ok = 0;
    eval { $doc->CalculateBonds; $ok = 1; };
    kv("CalculateBonds()", $ok ? "OK" : "FAILED");
    kv("UnitCell bonds (after CalculateBonds)", scalar(@{ $uc->Bonds }));
} else {
    kv("CalculateBonds()", "SKIPPED");
}

if ($REPORT_STYLE ne "tsv") {
    hr("=");
    print "[DONE]\n";
}

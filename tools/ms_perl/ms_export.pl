#!/usr/bin/perl
use strict;
use warnings;
use MaterialsScript qw(:all);
use Cwd qw(getcwd);
use File::Spec;

my $SUMMARY = "export_summary.txt";

sub list_subdirs {
    my ($dir) = @_;
    my @subs;
    opendir(my $dh, $dir) or return @subs;
    while (my $e = readdir($dh)) {
        next if $e eq "." || $e eq "..";
        my $p = File::Spec->catdir($dir, $e);
        push @subs, $p if -d $p;
    }
    closedir($dh);
    return @subs;
}

sub list_xsd_files_in_dir {
    my ($dir) = @_;
    my @files;
    opendir(my $dh, $dir) or return @files;
    while (my $e = readdir($dh)) {
        next unless $e =~ /\.xsd$/i;
        my $p = File::Spec->catfile($dir, $e);
        push @files, $p if -f $p;
    }
    closedir($dh);
    return @files;
}

# 尝试导出：按扩展名自动选择导出器；必要时可加格式名兜底
sub try_export_auto {
    my ($doc, $out) = @_;
    my $ok = 0;
    my $err = "";
    eval { $doc->Export($out); $ok = 1; };
    if (!$ok) { $err = $@; $err =~ s/[\r\n]+/ /g; }
    return ($ok, $err);
}

# ===== summary =====
my $sumdoc = Documents->New($SUMMARY);
$sumdoc->Append("XSD export summary (try CIF first; if fails export MOL)\n");
$sumdoc->Append("file\tfinal_format\tstatus\tmessage\toutput\n\n");

# ===== 找 xsd（cwd + 两层子目录）=====
my $cwd = getcwd();
my @candidates = ($cwd);
my @lvl1 = list_subdirs($cwd);
push @candidates, @lvl1;
foreach my $d1 (@lvl1) { push @candidates, list_subdirs($d1); }

my %seen;
my @all_xsd;
foreach my $dir (@candidates) {
    my @xsd = list_xsd_files_in_dir($dir);
    foreach my $p (@xsd) {
        next if $seen{$p}++;
        push @all_xsd, $p;
    }
}

if (!@all_xsd) {
    $sumdoc->Append("N/A\tN/A\tFAIL\tNo XSD found under Script cwd (up to 2 levels).\t\n");
    $sumdoc->Save();
    die "No XSD found.\n";
}

print "XSD count: " . scalar(@all_xsd) . "\n";

# ===== 主循环：先 CIF，失败再 MOL =====
foreach my $path (@all_xsd) {

    my ($vol, $dirs, $file) = File::Spec->splitpath($path);
    my $base = $file; $base =~ s/\.xsd$//i;

    my $status = "FAIL";
    my $fmt    = "";
    my $msg    = "";
    my $out    = "";

    eval {
        my $doc = $Documents{$file};
        if (!$doc) { $doc = Documents->Import($path); }
        $doc or die "Failed to get document: $file\n";

        # 1) 先尝试 CIF
        my $cif_out = "${base}.cif";
        my ($ok_cif, $err_cif) = try_export_auto($doc, $cif_out);

        if ($ok_cif) {
            $fmt = "CIF";
            $out = $cif_out;
            $status = "OK";
            $msg = "Exported as CIF";
        } else {
            # 2) CIF 失败 -> 导出 MOL
            my $mol_out = "${base}.mol";
            my ($ok_mol, $err_mol) = try_export_auto($doc, $mol_out);

            if ($ok_mol) {
                $fmt = "MOL";
                $out = $mol_out;
                $status = "OK";
                $msg = "CIF failed ($err_cif) ; exported as MOL";
            } else {
                $fmt = "N/A";
                $out = "N/A";
                die "CIF failed ($err_cif) ; MOL failed ($err_mol)";
            }
        }
    };

    if ($@) {
        $status = "FAIL";
        $msg = $@; $msg =~ s/[\r\n]+/ /g;
    }

    print "$file => $status ($fmt)\n";
    $sumdoc->Append("$file\t$fmt\t$status\t$msg\t$out\n");
}

$sumdoc->Save();
print "Done. Exported files + export_summary.txt are in Project root.\n";

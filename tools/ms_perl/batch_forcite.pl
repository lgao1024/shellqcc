#!/usr/bin/perl
use strict;
use warnings;
use MaterialsScript qw(:all);
use Cwd qw(getcwd);
use File::Spec;

# ====== Forcite 设置（按你的 Copy Script）======
my $FF       = "UFF4MOF";
my $QUALITY  = "Fine";
my $CHARGE   = "Use current";
my $P_EXT    = 0.0001;
my $CELL_opt  = "Yes";




# ====== 输出设置 ======
my $OutOptDir     = "Optimized_XSD";    # Project内输出目录（建议手动创建）
my $WRITE_SUMMARY = 1;                  # 不要summary就改成0
my $SummaryPath   = "$OutOptDir/summary.txt";  # summary 放进 Optimized_XSD

# ====== 工具函数：列子目录 ======
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

# ====== 工具函数：安全列目录内xsd（不用glob避免空格路径问题） ======
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

# ====== summary（可选） ======
my $sumdoc;
my $summary_actual_path = $SummaryPath;

if ($WRITE_SUMMARY) {
    # 先尝试写到 Optimized_XSD/summary.txt
    eval {
        $sumdoc = Documents->New($SummaryPath);
    };
    if ($@ || !$sumdoc) {
        # 如果 Optimized_XSD 不存在或不可写，就降级写到根目录
        $summary_actual_path = "summary.txt";
        $sumdoc = Documents->New($summary_actual_path);
    }

    $sumdoc->Append("Forcite batch geometry optimization summary\n");
    $sumdoc->Append("file\tstatus\tmessage\toptimized_xsd\n\n");
}

# ====== 自动搜索 XSD：script cwd + 两层子目录 ======
my $cwd = getcwd();
my @candidates = ($cwd);

my @lvl1 = list_subdirs($cwd);
push @candidates, @lvl1;

foreach my $d1 (@lvl1) {
    push @candidates, list_subdirs($d1);
}

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
    if ($WRITE_SUMMARY) {
        $sumdoc->Append("N/A\tFAIL\tNo XSD found under Script cwd (up to 2 levels).\t\n");
        $sumdoc->Save();
    }
    die "No XSD found.\n";
}

print "XSD count: " . scalar(@all_xsd) . "\n";

# ====== 主循环：只做优化 + 只保存最终opt.xsd到Optimized_XSD ======
foreach my $path (@all_xsd) {

    my ($vol, $dirs, $file) = File::Spec->splitpath($path);
    my $base = $file; $base =~ s/\.xsd$//i;

    my $status = "FAIL";
    my $msg    = "";
    my $opt_out = "";

    eval {
        # Project中已有同名文档则直接用，否则Import（避免重名Import失败）
        my $doc = $Documents{$file};
        if (!$doc) {
            $doc = Documents->Import($path);
        }
        $doc or die "Failed to get document for: $file\n";

        # Forcite 几何优化（你的Copy Script设置）
        Modules->Forcite->GeometryOptimization->Run($doc, Settings(
            Quality           => $QUALITY,
            CurrentForcefield => $FF,
            ChargeAssignment  => $CHARGE,
            ExternalPressure  => $P_EXT,
            OptimizeCell      => $CELL_opt, 
        ));

        # 保存最终优化结构到 Optimized_XSD（写不了就退回根目录）
        $opt_out = "$OutOptDir/${base}_opt.xsd";
        my $saved = 0;

        eval { $doc->SaveAs($opt_out); $saved = 1; };
        if (!$saved) {
            $opt_out = "${base}_opt.xsd";
            $doc->SaveAs($opt_out);
        }

        $status = "OK";
        $msg    = "Completed";
    };

    if ($@) {
        $status = "FAIL";
        $msg = $@; $msg =~ s/[\r\n]+/ /g;
    }

    print "$file => $status\n";
    if ($WRITE_SUMMARY) {
        $sumdoc->Append("$file\t$status\t$msg\t$opt_out\n");
    }
}

if ($WRITE_SUMMARY) {
    $sumdoc->Save();
    print "Summary saved to: $summary_actual_path\n";
}

print "Done. Optimized XSD saved to $OutOptDir (or root if folder not writable).\n";

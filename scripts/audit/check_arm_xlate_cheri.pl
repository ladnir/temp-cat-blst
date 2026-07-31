#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;
use File::Temp qw(tempdir);

my $repo_root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..', '..'));
my $translator = File::Spec->catfile($repo_root, 'blst', 'asm', 'arm-xlate.pl');
my $temp_root = tempdir('catblst-arm-xlate.XXXXXX', TMPDIR => 1, CLEANUP => 1);

my $input = <<'ASSEMBLY';
scvalue x10, x11, x12
scvalue x9, sp, x10
scbnds x13, x14, x15
scbnds x16, x17, #512
cadd x18, x19, x20
cadd x24, sp, #32
csub x21, x22, x23
ASSEMBLY

sub translate
{
    my ($flavour) = @_;
    my $output = File::Spec->catfile($temp_root, "$flavour.S");
    open(my $child, '|-', $^X, $translator, $flavour, $output)
        or die "can't start ARM translator for $flavour: $!\n";
    print {$child} $input or die "can't write ARM translator input: $!\n";
    close($child) or die "ARM translator failed for $flavour\n";

    open(my $generated, '<', $output)
        or die "can't read ARM translator output $output: $!\n";
    local $/;
    my $text = <$generated>;
    close($generated) or die "can't close ARM translator output $output: $!\n";
    $text =~ s/\r\n/\n/g;
    return $text;
}

sub require_lines
{
    my ($flavour, $text, @expected) = @_;
    for my $line (@expected) {
        die "$flavour translation is missing: $line\n"
            if $text !~ /^\Q$line\E$/m;
    }
}

my $cheri = translate('cheri64');
require_lines('cheri64', $cheri,
    "\tscvalue\tc10,c11,x12",
    "\tscvalue\tc9,csp,x10",
    "\tscbnds\tc13,c14,x15",
    "\tscbnds\tc16,c17,#512",
    "\tadd\tc18,c19,x20",
    "\tadd\tc24,csp,#32",
    "\tsub\tc21,c22,x23");

my $linux = translate('linux64');
require_lines('linux64', $linux,
    "\tadd\tx18,x19,x20",
    "\tsub\tx21,x22,x23");

print "ARM translator CHERI operand checks passed.\n";

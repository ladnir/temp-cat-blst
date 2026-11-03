#!/usr/bin/env perl

use strict;
use warnings;

my @file = <>;
my @structs;
my %seen;

sub require_rewrite
{
    my ($line_ref, $pattern, $replacement, $description) = @_;
    my $count = ($$line_ref =~ s/$pattern/$replacement/g);
    die "bindings trim expected exactly one $description rewrite, found $count\n"
        unless $count == 1;
}

for (my $i = 0; $i <= $#file; $i++) {
    $file[$i] =~ s/::std::os::raw::c_void/::core::ffi::c_void/g;
    $file[$i] =~ s/::std::os::raw::c_uint/u32/g;
    $file[$i] =~ s/::std::/::core::/g;
    # libclang joins adjacent Doxygen blocks with their comment delimiters.
    # Keep both blocks as one Rust doc string without exposing those markers.
    $file[$i] =~ s{\\n/\\n/\*\*}{\\n\\n}g;

    if ($file[$i] =~ m/pub\s+(?:struct|enum)\s+(\w+)/) {
        push @structs, $1;
    }

    next unless $file[$i] =~ m/^pub\s+struct\s+(\w+)/;
    my $name = $1;

    if ($name =~ m/^blst_p[12](?:_affine)?$/) {
        require_rewrite(\$file[$i - 1], qr/,\s*PartialEq/, "",
                        "$name PartialEq");
        $seen{point}++;
    } elsif ($name eq "blst_fp12") {
        require_rewrite(\$file[$i - 1], qr/,\s*Default/, "",
                        "$name Default");
        require_rewrite(\$file[$i - 1], qr/,\s*PartialEq/, "",
                        "$name PartialEq");
        $seen{fp12}++;
    } elsif ($name eq "blst_pairing_context" ||
             $name eq "blst_uniq_context" ||
             $name eq "blst_scalar_prepared_g2_batch") {
        for my $derive (qw(Copy Clone)) {
            require_rewrite(\$file[$i - 1], qr/,\s*$derive/, "",
                            "$name $derive");
        }
        $seen{opaque}++;
    } elsif ($name eq "blst_scalar") {
        require_rewrite(\$file[$i - 1], qr/,\s*Copy/, "",
                        "$name Copy");
        require_rewrite(\$file[$i - 1], qr/\)/, ", Zeroize)",
                        "$name Zeroize");
        splice @file, $i, 0, "#[zeroize(drop)]\n";
        $i++;
        $seen{scalar}++;
    }
}

my %expected = (point => 4, fp12 => 1, opaque => 3, scalar => 1);
for my $kind (sort keys %expected) {
    my $actual = $seen{$kind} // 0;
    die "bindings trim expected $expected{$kind} $kind types, found $actual\n"
        unless $actual == $expected{$kind};
}

print @file;

print << '___';
#[test]
fn bindgen_test_layout_blst_error() {
    assert_eq!(::core::mem::size_of::<BLST_ERROR>(), 4usize);
    assert_eq!(::core::mem::align_of::<BLST_ERROR>(), 4usize);
    assert_eq!(BLST_ERROR(0xffff_ffff).0, 0xffff_ffff);
}

#[test]
fn bindgen_test_layout_blst_uniq_result() {
    assert_eq!(::core::mem::size_of::<BLST_UNIQ_RESULT>(), 4usize);
    assert_eq!(::core::mem::align_of::<BLST_UNIQ_RESULT>(), 4usize);
    assert_eq!(BLST_UNIQ_RESULT(0xffff_ffff).0, 0xffff_ffff);
}

#[test]
fn bindgen_test_normal_types() {
    // from "Rust for Rustaceans" by Jon Gjengset
    fn is_normal<T: Sized + Send + Sync + Unpin>() {}
___
for (@structs) {
    print "    is_normal::<$_>();\n";
}
print "}\n";

close STDOUT;

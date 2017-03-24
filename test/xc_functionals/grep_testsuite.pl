#!/usr/bin/perl
use strict;
use warnings;
my %xc;
foreach my $filename (@ARGV) {
    next if ($filename =~ m#json|annotate# or not -f $filename);
    next unless open(my $fh, $filename);
    my $xc;
    my $exx;
    while (<$fh>) {
        chomp;
        if (m#^\s*Exchange-correlation\s*=\s*(.*)#) {
            print STDERR "xc already defined in $filename\n" if (defined $xc and ($xc ne $1));
            $xc = $1;
        } elsif (m#^\s*EXX-fraction\s*=\s*([0-9\.]+)\s*$#) {
            print STDERR "exx already defined in $filename\n" if (defined $exx and abs($exx-$1)>0.000001);
            $exx = $1;
        }
    }
    if (defined $xc) {
        my $xc_str = sprintf("%4.2f %s", 0, $xc);
        $xc{$xc_str}++;
        if (defined $exx and $exx > 0.000001) {
            my $xc_str = sprintf("%4.2f %s", $exx, $xc);
            $xc{$xc_str}++;
        }
    } else {
        print STDERR "no xc in $filename\n";
    }
}
my @sorti = sort { $xc{$a} <=> $xc{$b} } keys %xc;
foreach my $key (@sorti) {
    printf "%5d %s\n", $xc{$key}, $key;
}

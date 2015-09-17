#!/usr/bin/perl -w
use strict;
# perl HomecodingV1.0.pl blastFile fastaFile outputfulllengthcDNAFile outputextensionfullLengthcDNAFile 


use get_454codingv1_1;
my ($flcdna, $expect_flcdna)= get_454codingv1_1::codingext($ARGV[0], $ARGV[1]);
open FLCDNA, ">$ARGV[2]" || die "Can't open $ARGV[2]: $!";
open EFLCDNA, ">$ARGV[3]" || die "Can't open $ARGV[3]: $!";

print FLCDNA $flcdna;
print EFLCDNA $expect_flcdna;

close FLCDNA;
close EFLCDNA;


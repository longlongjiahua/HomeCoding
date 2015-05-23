#!/usr/bin/perl -w
use strict;
use get_454codingv1_1;
my ($flcdna, $expect_flcdna)= get_454codingv1_1::codingext($ARGV[0], $ARGV[1]);
open FLCDNA, ">$ARGV[2]";
open EFLCDNA, ">$ARGV[3]";

print FLCDNA $flcdna;
print EFLCDNA $expect_flcdna;

#!/use/bin/perl -w
use strict;
#May 2009 @ Yong Xu
#The algorithm was designed by the idea from 
#Salse et al  Plant Cell. 2008 Jan;20(1):11-24
##inter_species: Orthlogous pairs : 60% CIP, 70% CALP;
##IntraGenome Paralogous pairs :70% CIP, 70% CALP;


use cip_calp;
#$blastfile,$id_cutoff_inter,$cover_cutoff_inter,$id_cutoff_intra, $cover_cutoff_intra
my($intra, $inter)=cip_calp::cip_calp($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3], $ARGV[4]);

               cip_calp::orth_para($intra, $inter);

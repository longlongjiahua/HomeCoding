#!/usr/bin/perl -w
package get_blast_start_end;
use strict;
sub tiling{
      my @hsps=@_;  
      my $len_hsps=0;
      my @sorted=sort {$a->{hit_start}<=>$b->{hit_start}} @hsps;
      my @sorted_end=sort{$a->{hit_end}<=>$b->{hit_end}} @hsps;
      my $final_hit_end=$sorted_end[-1]->{hit_end};
      my $final_hit_start=$sorted[0]->{hit_start};
      #print "end- ", $final_end, "\n";
      #print join "\n", map {$_->{hit_start}." - ".$_->{hit_end}."-".$_->{similar}} @sorted; 
      #print "\n************************************************************************************\n";
      my @s_t=();
      my($five_ter, $three_ter)=get_segment(@sorted);
      #print $hitname, "\t", $hitlength, "\t","start  ", $five_ter, "\t", "end  ", $three_ter, "\n";
      my $len_seg=$three_ter-$five_ter+1;
      push @s_t, $five_ter."..".$three_ter;
      $len_hsps+=$len_seg;      
      while(1){
          if($three_ter<$final_hit_end){
                my @left=grep{$_->{hit_end}>$three_ter} @sorted;
               ($five_ter, $three_ter)=get_segment(@left);
	       my $len_seg=$three_ter-$five_ter+1;
	       push @s_t, $five_ter."..".$three_ter;
	       $len_hsps+=$len_seg;
              }
          else{last;}
         #print $hitname, "\t", $hitlength, "\t","start  ", $five_ter, "\t", "end  ", $three_ter, "\n";         
   } 
      my @sorted_query=sort {$a->{query_start}<=>$b->{query_start}} @hsps;
      my @sorted_query_end=sort{$a->{query_end}<=>$b->{query_end}} @hsps;
      my $final_query_start=$sorted_query[0]->{query_start};
      my $final_query_end=$sorted_query_end[-1]->{query_end};
  return $final_hit_start, $final_hit_end, $len_hsps,$final_query_start, $final_query_end;        
 }# main sub;

sub get_segment{

      my @sorted=@_;
      my $key1=shift @sorted; 
      my $five_ter=$key1->{hit_start};
      my $three_ter=$key1->{hit_end};
      foreach my $key (@sorted){ 
            my $start=$key->{hit_start};
	    my $end=$key->{hit_end};
            #print $start, "\t", $end, "\n";
	    $five_ter=$start if($start<$five_ter);
	    $three_ter=$end if($three_ter<$end&&$start<$three_ter); 
        }
   return $five_ter, $three_ter;   
 }
1;

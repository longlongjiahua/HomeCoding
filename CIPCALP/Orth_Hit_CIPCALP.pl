#!usr/bin/perl -w
use strict;

#input blastafile identical_cutoff
##inter_species: Orthlogous pairs : 60% CIP, 70% CALP;
##IntraGenome Paralogous pairs :70% CIP, 70% CALP;
#calp ->coverage;
use warnings; 
my %inter;
my %intra;
my($blastfile, $out, $id_cutoff_inter,$cover_cutoff_inter)=@ARGV;
open OUT, ">$out"; my %uni=();
my @cipcalp=();##change
use Bio::SearchIO; 
my $blast= new Bio::SearchIO (-file=>$blastfile, -format=>'blast');
while (my $result = $blast->next_result){ 
       my $queryname=$result->query_accession;
       my $querylength=$result->query_length;
       my $jud_orth=0;
       #print $queryname, "\n";	
       my $good_hits=0;       
       while (my $hit = $result->next_hit){#######
	     my $hitlength=$hit->length;
	     my $hitname=$hit->accession;
	     my @query_inds=(); my @hit_inds=(); 
	     my @t_length=();
	     my $i=0;
	     next if($queryname eq $hitname);
	     while (my $hsp = $hit->next_hsp){	
	         $i++;
		 my $percent_id=$hsp->percent_identity;
		 last if($percent_id<70&&$i==1);			 					
	         my $hsplength=$hsp->length('query');	      
	         my $similar=$hsp->percent_identity;
		 my $hit_strand=$hsp->strand('hit');
		 my $hit_start=$hsp->start('hit');
		 my $hit_end=$hsp->end('hit');
		 my $query_start=$hsp->start('query');
		 my $query_end=$hsp->end('query');
		 push @query_inds,$hsp->seq_inds('query','identical');
		 push @hit_inds,$hsp->seq_inds('hit','identical');
		 #my $cover=$hsplength/$querylength;
		 #print $query_start,"\t",$query_end,"\n";		 
		 push @t_length, ($hit_start..$hit_end);		 
	     }
	 my $jud=   judge(\@hit_inds, \@t_length, $hitname,$queryname,  $hitlength,
                   $id_cutoff_inter,$cover_cutoff_inter);
	  if($jud==1){$good_hits++; print $queryname, "\t", $hitname, "\n";}
         $uni{$hitname}++ if($jud==1);	   
         $jud_orth=1  if($jud==1);
     }
     print OUT $queryname."\n" if($jud_orth==0);
    # print $queryname."\t", $good_hits, "\n";
     
}
#while(my ($k, $v)=each %uni){
my $num=keys %uni;
print "Total hit: ", $num, "\n";
 


 
sub judge {
     my ($query_ind, $t_len, $hitname, $queryname, $querylength,
                   $id_cutoff_inter,$cover_cutoff_inter)=@_;
            my @query_inds = @{$query_ind};
	    my @t_length   = @{$t_len};
	    my %seen1=(); 
	    my %seen2=();	    
	    my @sum_ids=grep{!$seen1{$_}++} @query_inds;
	    my @als=grep{!$seen2{$_}++} @t_length;
	    if($#sum_ids<0||$#als<0){
	        #print $queryname, "\t", $hitname, "\n";
		return 0;
	     }  
	    #print join("t", @als), "\n";
	    my $t_ids=1+$#sum_ids;	
	    my $t_al=1+$#als; 	    
            my $cip=$t_ids/$t_al*100;	    
            my $calp=$t_al/$querylength*100;
	    #print $query_sim, "\t", $hit_sim,"\t",$hitname, "\n";
            return 0 if($cip<$id_cutoff_inter);
            return 0 if($calp<$cover_cutoff_inter);
	    #print $queryname."\t".$hitname."\t".$cip."\t".$calp, "\n";
            return 1;
}    
 

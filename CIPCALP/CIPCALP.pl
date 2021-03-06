#!usr/bin/perl -w
use strict;

#input blastafile identical_cutoff
##inter_species: Orthlogous pairs : 60% CIP, 70% CALP;
##IntraGenome Paralogous pairs :70% CIP, 70% CALP;
#calp ->coverage;
use warnings; 
my %inter;
my %intra;
my($blastfile, $out, $id_cutoff_inter,$cover_cutoff_inter,$id_cutoff_intra, $cover_cutoff_intra)=@ARGV;
open OUT, ">$out";
my @cipcalp=();##change
use Bio::SearchIO;
my $blast= new Bio::SearchIO (-file=>$blastfile, -format=>'blast');
while (my $result = $blast->next_result){ 
       my $queryname=$result->query_accession;
       my $querylength=$result->query_length;
       my $jud_orth=0;
       #print $queryname, "\n";	       
       while (my $hit = $result->next_hit){#######
	     my $hitlength=$hit->length;
	     my $hitname=$hit->accession;
	     my @query_inds=();
	     my @t_length=();
	     my $i=0;
	     next if($queryname eq $hitname);
	     while (my $hsp = $hit->next_hsp){	
	         $i++;
		 my $percent_id=$hsp->percent_identity;
		 last if($percent_id<60&&$i==1);			 					
	         my $hsplength=$hsp->length('query');	      
	         my $similar=$hsp->percent_identity;
		 my $hit_strand=$hsp->strand('hit');
		 my $hit_start=$hsp->start('hit');
		 my $hit_end=$hsp->end('hit');
		 my $query_start=$hsp->start('query');
		 my $query_end=$hsp->end('query');
		 push @query_inds,$hsp->seq_inds('query','identical');
		 #my $cover=$hsplength/$querylength;
		 #print $query_start,"\t",$query_end,"\n";		 
		 push @t_length, ($query_start..$query_end);		 
	     }
	 my $jud=   judge(\@query_inds, \@t_length, $hitname, $queryname, $querylength,
                   $id_cutoff_inter,$cover_cutoff_inter,$id_cutoff_intra, $cover_cutoff_intra);
         $jud_orth=1 if($jud==1);
     }
     print OUT $queryname."\n" if($jud_orth==0);
} 
sub judge {
     my ($query_ind, $t_len, $hitname, $queryname, $querylength,
                   $id_cutoff_inter,$cover_cutoff_inter,$id_cutoff_intra, $cover_cutoff_intra)=@_;
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
	    my $query_sim=substr($queryname, 0,2);
	    my $hit_sim=substr($hitname, 0, 2);
	    #print $query_sim, "\t", $hit_sim,"\t",$hitname, "\n";
	    if($query_sim eq $hit_sim){
                return 0 if($cip<$id_cutoff_intra);
                return 0 if($calp<$cover_cutoff_intra);
	        print $queryname."\t".$hitname."\t".$cip."\t".$calp, "\n";
		return 1;
                #push @cipcalp, {querynam=>$queryname, hitname=>$hitname, cip=>$cip, calp=>$calp};
		#push @{$intra{$queryname}},$hitname;
	     }
	    if($query_sim ne $hit_sim){
                #return if($cip<$id_cutoff_inter);
                #return if($calp<$cover_cutoff_inter);
	        #print $queryname."\t".$hitname."\t".$cip."\t".$calp, "\n";
                #push @cipcalp, {querynam=>$queryname, hitname=>$hitname, cip=>$cip, calp=>$calp};
		#push @{$inter{$queryname}},$hitname;
	     }
 }

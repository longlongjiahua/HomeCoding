package cip_calp;
#input blastafile identical_cutoff
##inter_species: Orthlogous pairs : 60% CIP, 70% CALP;
##IntraGenome Paralogous pairs :70% CIP, 70% CALP;
use strict;
#calp ->coverage;
use warnings; 
my %inter;
my %intra;
sub cip_calp{
  my($blastfile,$id_cutoff_inter,$cover_cutoff_inter,$id_cutoff_intra, $cover_cutoff_intra)=@_;
  my @cipcalp=();##change
  use Bio::SearchIO;
  my $blast= new Bio::SearchIO (-file=>"$blastfile", -format=>'blast');
  while (my $result = $blast->next_result){ 
       my $queryname=$result->query_accession;
       my $querylength=$result->query_length;
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
	    my %seen1=(); 
	    my %seen2=();	    
	    my @sum_ids=grep{!$seen1{$_}++} @query_inds;
	    	    
	    my @als=grep{!$seen2{$_}++} @t_length;
	    if($#sum_ids<0||$#als<0){
	        #print $queryname, "\t", $hitname, "\n";
		next;
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
                next if($cip<$id_cutoff_intra);
                next if($calp<$cover_cutoff_intra);
	        #print $queryname."\t".$hitname."\t".$cip."\t".$calp, "\n";
                #push @cipcalp, {querynam=>$queryname, hitname=>$hitname, cip=>$cip, calp=>$calp};
		push @{$intra{$queryname}},$hitname;
	     }
	    if($query_sim ne $hit_sim){
                next if($cip<$id_cutoff_inter);
                next if($calp<$cover_cutoff_inter);
	        #print $queryname."\t".$hitname."\t".$cip."\t".$calp, "\n";
                #push @cipcalp, {querynam=>$queryname, hitname=>$hitname, cip=>$cip, calp=>$calp};
		push @{$inter{$queryname}},$hitname;
	     }
     }
   
  }
 
  #return @cipcalp;
  return \%intra, \%inter; 
    
} 
  sub orth_para{
       my ($intra_h, $inter_h)=@_;
       open OS_SOLO, ">os_single";
       open BRADI_SOLO, ">bradi_single";
       open OS_ORTH, ">os_orth";
       open BRADI_ORTH, ">bradi_orth";
       
       my $os_single=0;
       my $bradi_single=0;
       my $os_orth=0;
       my $bradi_orth=0;              
       my %intra=%{$intra_h};
       my %inter=%{$inter_h};       
       foreach my $name (keys %intra){
           # print $name, ${$intra{$name}}[1], "\n";
              next if($intra{$name}[1]);
	      #print $name, ${$intra{$name}}[0], "\n";
	      if($name=~/Os/){
		     $os_single++;
                     print OS_SOLO  ($name,"\t",$intra{$name}[0],"\n");
                  }   
                if($name=~/Bradi/){
		      $bradi_single++;
		      print BRADI_SOLO ($name,"\t",$intra{$name}[0],"\n");
                  }
           
	}  
       foreach my $name1 (keys %inter){
            #print $name1, @{$inter{$name1}}, "\n";
            next if($inter{$name1}[1]);
	        if($name1=~/Os/){
		     $os_orth++;
                     print OS_ORTH ($name1,"\t",$inter{$name1}[0],"\n");
                  }   
                if($name1=~/Bradi/){
		      $bradi_orth++;
		       print BRADI_ORTH ($name1,"\t",$inter{$name1}[0],"\n"); 
                  }
          
	}     
	   
 }	   
1; 

#plus strand ......A(7)TG-----TGA(21)........;
#minusstrand ......T(7)CA----CAT(21).........;
#perl get_partial.pl /home/xuyong1/IDASM/Divide_Stitch_merge/V5.2/V52P3ToGrassFL.blastn /home/xuyong1/IDASM/Divide_Stitch_merge/V5.2/V52P3_NoFLCorrect_16250.fsa >/home/xuyong1/IDASM/Divide_Stitch_merge/V5.2/partial

package for_partial;
###this packgae get the blastfile and fastafile, predict 
use strict;
use Bio::SearchIO;
use warnings; 
use get_blast_start_end;

sub codingext { 
  my($blastfile, $fastafile)=@_;
  my %name_seq=name_seq($fastafile);
  my $blast= new Bio::SearchIO (-file=>$blastfile, -format=>'blast');
  my $kk=0;
  while (my $result = $blast->next_result){   
       my $queryname=$result->query_accession;       
       my $querylength=$result->query_length;
       my $m=0;
       my $n=0;
       my @hsps=();
       my $ii=0;      	       
       while(my $hit=$result->next_hit){#######
             $ii++;					
	     last if($ii>1);
	     my $hitlength=$hit->length;	     	     
	     my $hitname=$hit->accession;
	     my $seq=$name_seq{$queryname};
	     my $strand=0;
	     my $i_strand=0;
	     while (my $hsp = $hit->next_hsp){
	         $i_strand++;					 					
	         my $hsplength=$hsp->length('query');	       
	         my $similar=$hsp->percent_identity;
		 next if($similar<40);
		 my $hit_strand=$hsp->strand('hit');
		 my $hit_start=$hsp->start('hit');
		 my $hit_end=$hsp->end('hit');
		 my $query_start=$hsp->start('query');
		 my $query_end=$hsp->end('query');		 
		 my $cover=$hsplength/$hitlength;
		 $strand=$hsp->strand('query') if($i_strand==1);
		 push @hsps, {hit_start=>$hit_start, hit_end=>$hit_end, query_start=>$query_start, 
		  query_end=>$query_end, similar=>$similar};
	      }#parsing one hit finished
	     last if($#hsps<0);
	     my($fhit_start, $fhit_end,$len_hsps, $fquery_start, $fquery_end)=get_blast_start_end::tiling(@hsps);
#########	     #####if($len_hsps/$hitlength>=0.5){
             if($len_hsps/$hitlength>=0.5||$len_hsps>50){
	              my ($head, $fas)=xts($fhit_start, $fhit_end, $hitlength, 
	                                 $fquery_start, $fquery_end, $querylength,$strand, $seq, $queryname);
		      if($head){
		            $head.=" hitname:".$hitname." HitLength:". $hitlength. " HitHSP:".$fhit_start."_".$fhit_end."\n";
			    #strand:1 queryLength:929 expect:115_360 QueryHSP:115_360 extension:115_360 hitname:AK250501 HitHSP:
	
			    print ">".$head, $fas, "\n";
		      }	 
		}
     }
 }#parse blast result; 
} #end of main sub;
sub xts{
   my($hit_start, $hit_end, $hitlength, $query_start, $query_end, $querylength,$strand,$seq,$queryname)=@_;
                 my $five_ter=0;
		 my $three_ter=0;
		 my $expect_start=0;
		 my $expect_end=0;
		 my $fasseq="";
		 my $five_expect=3*($hit_start-1);#absolute length
		 my $three_expect=3*($hitlength-$hit_end);#hit seq as reference
		 if($strand==-1){#minus strand
		    $expect_start=$query_start-$three_expect;#absolute length
		    $expect_end=$query_end+$five_expect;		    
		    if($query_start>$three_expect){###here utr respresent terminal
		         my $expect_pos=$query_start-$three_expect;
	                 $five_ter=extension_reverse($hit_start, $hit_end, $query_start, $query_end,
		                                          $seq, $querylength, $strand,$expect_pos, $queryname);
	             }
		    if($query_start<=$three_expect){
		        my $uncom=$three_expect-$query_start;
		       $five_ter="uncomplete in 5' terminal ".$uncom." nt ";
		     }
		    if ($hit_start>1&&($querylength-$query_end)>$five_expect){
		         my $expect_pos=$query_end+$five_expect;
	                 $three_ter=extension_forward($hit_start, $hit_end, $query_start, $query_end,$seq,  
		                                                      $querylength,$strand,$expect_pos, $queryname);
	               } 
		    if(($querylength-$query_end)<=$five_expect){ 
		        my $uncom=$five_expect-($querylength-$query_end);
			$three_ter="uncomplete in 3' terminal ".$uncom." nt ";
		       }
		}####strand minus; 
		if ($strand==1){		    
		    $expect_start=$query_start-$five_expect;
		    $expect_end=$query_end+$three_expect;
		     if($hit_start>1&&($query_start>$five_expect)){
		         my $expect_pos=$query_start-$five_expect;
	                 $five_ter=extension_reverse($hit_start, $hit_end, $query_start, $query_end,
		                                          $seq, $querylength, $strand,$expect_pos, $queryname);
	                }
		     if($hit_start>1&&($query_start<=$five_expect)){
		         my $uncom=abs($expect_start);
		         $five_ter="uncomplete in 5' terminal ".$uncom."nt ";
		       }
		     if(($querylength-$query_end)>$three_expect){
		         my $expect_pos=$query_end+$three_expect;
	                 $three_ter=extension_forward($hit_start, $hit_end, $query_start, $query_end,$seq,  
		                                                      $querylength, $strand,$expect_pos, $queryname);
	               }
		     if(($querylength-$query_end)<$three_expect){  
		         my $uncom=$expect_end-$query_end;
		         $three_ter="uncomplete in 3' terminal ".$uncom."nt ";
		      }
	       } ###plus strand;
	     if($hit_start==1){
                  $five_ter=$query_start if($strand==1);
                  $three_ter=$query_end if($strand==-1);
              }
	     if($hit_end==$hitlength){
                   $five_ter=$query_start if($strand==-1);
                   $three_ter=$query_end if($strand==1);
              }
	     if($five_ter=~/uncomplete/||$three_ter=~/uncomplete/){
	        my $flcds=$expect_end-$expect_start+1;
	        my $five_act=0;
	        my $three_act=0;	     
	        $five_act=($five_ter=~/uncomplete/ ? 1 :$expect_start);
	        $three_act=($three_ter=~/terminal/ ? $querylength : $expect_end);
	         my $length_act=$three_act-$five_act+1;
	        my $ratio=$length_act/$flcds * 100;
	        $fasseq= ">".$queryname." ".$strand. " QueryLength:". $querylength ." expect:".$expect_start."_".$expect_end." "."HSP:".$query_start
	         ."_".$query_end." "."extension: ".$five_ter."||".$three_ter."  ratio:".$ratio."%"."\n".$seq."\n";
	       my $head = $queryname." strand:". $strand." QueryLength:". $querylength ." expect:". $expect_start. "_". $expect_end." QueryHSP:". $query_start
	         ."_".$query_end. " extension:".$five_ter."||".$three_ter."  ratio:".$ratio."%";
	       return ($head, $seq);
	    }
}
#########################
sub name_seq{ 
 use Bio::SeqIO;
 my %name_seq=();
 my $infasta= new Bio::SeqIO (-file=>$_[0], -format=>"fasta");
 while(my $fasta=$infasta->next_seq){
    my $name=$fasta->display_id;
    my $seq=$fasta->seq();   
   $name_seq{$name}=$seq;
  } 
 return (%name_seq);
}  
########################################## 
sub extension_forward{  
 my ($hit_start, $hit_end, $query_start, $query_end, $seq, $querylength,$strand,$expect_pos, $queryname)=@_;    
 my %pos=();
 my $scan_start=$expect_pos-30;
 my $scan_end=length($seq);
 $scan_end=$expect_pos+30 if(($expect_pos+30)<=length($seq)); 
 if($strand==1){
     for(my $i=$scan_start;$i<=$scan_end;$i+=3){
         my $tri_nt=substr($seq, $i,3);
	#print $i, "\t",$tri_nt, "\n";
	if($tri_nt=~/TAG|TGA|TAA/){
           $pos{$i}=abs($expect_pos-$i);	   
	}
    }
    if(keys %pos){    
        my @sorted= sort {$pos{$a} <=> $pos{$b}} (keys %pos);
        return ($sorted[0]+3);
    }    
 }  
 if($strand==-1){
    my %pos=();    
    my $scan_start=$expect_pos-30;
    my $scan_end=length($seq);
    $scan_end=$expect_pos+30 if(($expect_pos+30)<=length($seq));            
    for(my $i=$scan_start;$i<=$scan_end;$i+=3){
        my $tri_nt=substr($seq, $i,3);
	if($tri_nt=~/CAT/){
           $pos{$i}=abs($expect_pos-$i);	   
	}
    }
    if(keys %pos){    
        my @sorted= sort {$pos{$a} <=> $pos{$b}} (keys %pos);
        return ($sorted[0]+3);
    }    
  }
 } 
######################
sub extension_reverse {
    my ($hit_start, $hit_end, $query_start, $query_end,$seq,$querylength, $strand,$expect_pos, $queryname)=@_;    
    my %pos=();          
    my $int=int($expect_pos/3);
    my $scan_start=$expect_pos-(3*$int);
    $scan_start=$expect_pos-30 if(($expect_pos-30)>=0);
 if($strand==1){ 
     for(my $i=$scan_start;$i<=($expect_pos+30);$i+=3){
         my $tri_nt=substr($seq, $i-1,3);
	 if($tri_nt=~/ATG/){
           $pos{$i}=abs($expect_pos-$i);
	 }
     }
     #print "start     ",$expect_pos, "\n", join("\t", keys %pos),"\n";
    if(keys %pos){
       my @sorted= sort {$pos{$a} <=> $pos{$b}} (keys %pos);
       return ($sorted[0]);
    }
 } 
if($strand==-1){     
    for(my $i=$scan_start;$i<=($expect_pos+30);$i+=3){
        my $tri_nt=substr($seq, $i-1,3);
	if($tri_nt=~/CTA|TCA|TTA/){
           $pos{$i}=abs($expect_pos-$i);
	}
    }
    if(keys %pos){
       my @sorted= sort {$pos{$a} <=> $pos{$b}} (keys %pos);
       return ($sorted[0]);
    } 	       
 }
}
################
sub rev_comp {
     my $seq=$_[0];
     my $len=length($seq);
     my $five_ter=$_[1];
     my $three_ter=$_[2];
     my $r_seq=reverse($seq);
     
     $r_seq=~tr/ATGC/TACG/;
     
     my $start=$len-$three_ter+1;
     my $end=$len-$five_ter+1;
     return ($r_seq, $start, $end);    
}

1;

package get_454codingv1_1; 
#new version July 2010
###this packgae get the blastfile and fastafile, predict 
#added May 7th: looking for the second if the first is not complete ORF;
#added May 7th: if the part of query seqeuence was used two times, seleted the highest identical one
use strict;
#plus strand ......A(7)TG-----TGA(21)........;
#minusstrand ......T(7)CA----CAT(21).........;
use Bio::SearchIO;
use warnings; 
use get_blast_start_end;
sub codingext { 
  my $ffna="";
  my $expect_fna="";
  my($blastfile, $fastafile)=@_;
  my %name_seq=name_seq($fastafile);
  my $blast= new Bio::SearchIO (-file=>"$blastfile", -format=>'blast');
  my $ii=0;
  while (my $result = $blast->next_result){  
       #$ii++;last if($ii>20); 
       my $queryname=$result->query_accession;       
       my $querylength=$result->query_length;
       #next unless($queryname eq "Contig2");       
       print $queryname,"\t", $querylength, "\n";
       my $m=0;
       my $n=0;
       my @hsps=();
       #print $queryname, "\n"; 
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
		 next if($similar<30);
		 my $hit_strand=$hsp->strand('hit');
		 my $hit_start=$hsp->start('hit');
		 my $hit_end=$hsp->end('hit');
		 my $query_start=$hsp->start('query');
		 my $query_end=$hsp->end('query');		 
		 my $cover=$hsplength/$hitlength;
		 $strand=$hsp->strand('query') if($i_strand==1);
                 my $a_strand=$strand=$hsp->strand('query');
                 next if($a_strand ne $strand);  ###added on May 7th, to get rid of the partial with reverse direction
		 push @hsps, {hit_start=>$hit_start, hit_end=>$hit_end, query_start=>$query_start, 
		              query_end=>$query_end, similar=>$similar};
		}#parsing one hit finished
		last if($#hsps<0);
	        my($fhit_start, $fhit_end,$len_hsps, $fquery_start, $fquery_end)=get_blast_start_end::tiling(@hsps);
	        if($len_hsps/$hitlength>=0.5){
	        #print "tiling output:  ", "\t","strand: ", $strand,  "\t hit_end: ",$fhit_end, "\t hitlength: ", $hitlength,"\n";
	        my($tt,$fna)=xts($fhit_start, $fhit_end, $hitlength, 
	                         $fquery_start, $fquery_end, $querylength,$strand, $seq, $queryname, $hitname);		
		if($tt eq "true"){
	               #print $fna, "\n";			
	               $ffna.=$fna;
		}
		else {
		        $expect_fna.=$fna;
		}
	       }	           	                         
     }
 }#parse blast result; 
 return ($ffna, $expect_fna);
} #end of main sub;

sub xts{
  my($hit_start,$hit_end, $hitlength, $query_start, $query_end, $querylength,$strand,$seq,$queryname,$hitname)=@_;
  #print "input xts: ", join("\t", ($hit_start, $hit_end, $hitlength, $query_start, $query_end, $querylength,$strand)), "\n";
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
		    if ($query_start>$three_expect){###here utr respresent terminal
		         my $expect_pos=$query_start-$three_expect;
	                 $five_ter=extension_reverse($hit_start, $hit_end, $query_start, $query_end,
		                                          $seq, $querylength, $strand,$expect_pos, $queryname);
	                } 
		    if ($hit_start>1&&($querylength-$query_end)>$five_expect){
		         my $expect_pos=$query_end+$five_expect;
	                  $three_ter=extension_forward($hit_start, $hit_end, $query_start, $query_end,$seq,  
		                                                      $querylength,$strand,$expect_pos, $queryname);
                         #print $queryname, $strand, "\t",$hitname,"\t",$hit_end, "\t", $three_ter, "\n";
	               } 
		}####strand minus; 
		if ($strand==1){		    
		    $expect_start=$query_start-$five_expect;
		    $expect_end=$query_end+$three_expect;
		     if ($hit_start>1&&($query_start>$five_expect)){
		         my $expect_pos=$query_start-$five_expect;
	                 $five_ter=extension_reverse($hit_start, $hit_end, $query_start, $query_end,
		                                          $seq, $querylength, $strand,$expect_pos, $queryname);
                         #print $queryname, "\t",$hitname,"\t",$hit_start, "\t", $five_ter, "\n";
	                } 
		    if (($querylength-$query_end)>$three_expect){
		         my $expect_pos=$query_end+$three_expect;
		         
	                 $three_ter=extension_forward($hit_start, $hit_end, $query_start, $query_end,$seq,  
		                                                      $querylength, $strand,$expect_pos, $queryname);
                         #print $queryname, "\t",$hitname,"\t",$hit_end, "\t", $three_ter, "\n";
	               } 
		} ###plus strand;
               #print $queryname," ",$five_ter,"\t", $three_ter, "\n";
	      if($hit_start==1){
                  $five_ter=$query_start if($strand==1);
                  $three_ter=$query_end if($strand==-1);
              }
	      if($hit_end==$hitlength){
                      $five_ter=$query_start if($strand==-1);
                      $three_ter=$query_end if($strand==1);
               }

             if($five_ter>0&&$three_ter>0){  	        
                $fasseq= ">". $queryname."\t". "strand:" .$strand. " QueryLength:" .$querylength . " expect:".$expect_start. "_". 
		$expect_end. " QueryHSP:" .$query_start ."_". $query_end." extension:". $five_ter. "_". $three_ter. " hitname:". 
		$hitname." HitLength:".$hitlength. " HitHSP:".$hit_start."_".$hit_end."\n".$seq."\n";
		#print $fasseq;
	      return ("true",$fasseq);      
              }
	     else{
	            if($expect_start>0&&$expect_end<$querylength){
                $fasseq= ">". $queryname."\t". "strand:" .$strand. " QueryLength:" .$querylength . " expect:".$expect_start. "_". 
		$expect_end. " QueryHSP:" .$query_start ."_". $query_end." extension:". $five_ter. "_". $three_ter. " hitname:". 
		$hitname." HitLength:".$hitlength. " HitHSP:".$hit_start."_".$hit_end."\n".$seq."\n";
		#print $fasseq;

	                }	       	        
                return ("expect",$fasseq);
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
  #print ">", $name, "\n",$finalseq,"\n"; 
  } 
 return (%name_seq);
}  
########################################## 
sub extension_forward{  
 my ($hit_start, $hit_end, $query_start, $query_end, $seq, $querylength,$strand,$expect_pos, $queryname)=@_;    
 my %pos=();
 my $range =45;
 my $scan_start=$expect_pos-$range;
 my $scan_end=length($seq);
 $scan_end=$expect_pos+$range if(($expect_pos+$range)<=length($seq)); 
 if($strand==1){
     for(my $i=$scan_start;$i<=$scan_end;$i+=3){
        #for(my $i=119;$i<=200;$i+=3){
        my $tri_nt=substr($seq, $i,3);
	#print $i, "\t",$tri_nt, "\n";
	if($tri_nt=~/TAG|TGA|TAA/){
           $pos{$i}=abs($expect_pos-$i);
           last;	   
	}
    }
    if(keys %pos){    
        my @sorted= sort {$pos{$a} <=> $pos{$b}} (keys %pos);
        #print "start     ",$expect_pos, "\n", join("\t", keys %pos), "\n",join("\t", @sorted), "end   \n";
        return ($sorted[0]);
    }    
 }  
 if($strand==-1){
    my %pos=();    
    my $scan_start=$expect_pos-$range;
    my $scan_end=length($seq);
    $scan_end=$expect_pos+$range if(($expect_pos+$range)<=length($seq));            
    for(my $i=$scan_start;$i<=$scan_end;$i+=3){
        my $tri_nt=substr($seq, $i,3);
	#print $tri_nt, "\n";
	if($tri_nt=~/CAT/){
           $pos{$i}=abs($expect_pos-$i);	   
	}
    }
    if(keys %pos){    
        my @sorted= sort {$pos{$a} <=> $pos{$b}} (keys %pos);
        #print "start     ",$expect_pos, "\n", join("\t", keys %pos), "\n",join("\t", @sorted), "end   \n";
        return ($sorted[0]+3);
    }    
  }
 } 
######################
sub extension_reverse {
    my ($hit_start, $hit_end, $query_start, $query_end,$seq,$querylength, $strand,$expect_pos, $queryname)=@_;    
    my %pos=();
    my $range =39;          
    my $int=int($expect_pos/3);
    my $scan_start=$expect_pos-(3*$int);
    $scan_start=$expect_pos-$range if(($expect_pos-$range)>=0);
 if($strand==1){ 
     for(my $i=($expect_pos+$range);$i>=$scan_start;$i-=3){
         my $tri_nt=substr($seq, $i-1,3);
	 if($tri_nt=~/ATG/){
           $pos{$i}=abs($expect_pos-$i);
	 }
     }
     #print "start     ",$expect_pos, "\n", join("\t", keys %pos),"\n";
    if(keys %pos){
       my @sorted= sort {$pos{$a} <=> $pos{$b}} (keys %pos);
       print "start     ",$expect_pos, "\n", join("\t", keys %pos), "\n",join("\t", @sorted), "end   \n";
       return ($sorted[0]);
    }
 } 
if($strand==-1){     
    for(my $i=$scan_start;$i<=($expect_pos+$range);$i+=3){
        my $tri_nt=substr($seq, $i-1,3);
	if($tri_nt=~/CTA|TCA|TTA/){
           $pos{$i}=abs($expect_pos-$i);
	}
    }
    #print "start     ",$expect_pos, "\n", join("\t", keys %pos),"\n";
    if(keys %pos){
       my @sorted= sort {$pos{$a} <=> $pos{$b}} (keys %pos);
       #print "start     ",$expect_pos, "\n", join("\t", keys %pos), "\n",join("\t", @sorted), "end   \n";
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

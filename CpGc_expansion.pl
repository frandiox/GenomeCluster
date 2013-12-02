#!/usr/bin/perl -w


###########################################################################
###########################################################################

#              ********************************
#              **  CpGcluster (version 2.0)  **
#			         **		      stand-alone		     **
#              ********************************

#   Laboratorio de Genómica Evolutiva y BioInformática
#    Universidad de Granada, Departamento de Genética

#              Web: http://bioinfo2.ugr.es
#              CGI: http://bioinfo2.ugr.es/CpGcluster


#  For questions, feedback, etc please contact to: José L. Oliver (oliver@ugr.es)
#                                                  Michael Hackenberg (mlhack@gmail.com)
#                                                  Francisco Dios (frankdiox@gmail.com)

# To see the options of the algorithm, please launch CpGcluster without any command line arguments
# If you use CpGcluster, please cite...

############################################################################
############################################################################
use strict;


my ($seqst, $ID, $seqlen, $seqlenbruto, $elemnr, $prob, $output, $assembly_dir,
  $d, $plimit, $getd, @dd, @getd_hash ,@dist_n,$chrom_intersec,$genome_intersec,
  $bed_file, $nBedFile, %n_coords, %chrom_bed, %all_dist_n, %n_chrom_index, $elemmn);

#################################
### Parameters ##################
#################################
my $maxN = 0; ## maximal number of Ns
#################################

&GetDefault();


#my @f_assembly = split(/\//,$assembly_dir);
#my @f1_assembly = split(/\./,$f_assembly[$#f_assembly]);


my @dist_all; # holds distances for the genome
my $elemnr_genome; # number of genomic elements in genome
my $length_genome; # number of dinucleotides in genome
my $elmean_genome; # genomic element mean length in genome


####


print "      ***             Getting    Coordinates            ***\n";
GetBED();


#die "tiriri";

####

#if(opendir (my $DIR,$assembly_dir)){
	$elemnr_genome = 0;
	$length_genome = 0;
  $elmean_genome = 0;

	my $output_aux = $output;
	print "\n      --- Single chromosome:\n";

	for my $chrom (keys %chrom_bed){ ### MAIN LOOP

		print "\n      ".$chrom."\n";

		$output_aux = $output;


		#my @cod = &GetCoords($seqst,"CG"); #coordenadas de inicio de cada elemento (solapa con la primera posición del elemento)
		#my @cod = GetBED($assembly_dir."/".$chrom.".bed");

    #my $cod = $chrom_bed{$chrom}{'cStart'}; # $cod es un ARRAY REF


#open (S,">seqst.txt") or die "Can't open pinchado.txt";
#open (ID,">id.txt") or die "Can't open pinchado.txt";
#print S $seqst;
#print ID $ID;


#open (P,">pinchado.txt") or die "Can't open pinchado.txt";
#foreach my $elem (@cod){
	#print P $elem;


    print "      ***             Calculating   Seq   Features           ***\n";
    my ($num_dinuc,$last_chunk);
    $elemnr = $#{ $chrom_bed{$chrom}{'cStart'} } + 1; #numero de entradas almacenadas para este chrom
    my $sum1 = 0;
    for (my $i = 0; $i < $elemnr; $i++){
      $sum1 += $chrom_bed{$chrom}{'cEnd'}[$i] - $chrom_bed{$chrom}{'cStart'}[$i];
    }
    my $elmean = int(($sum1/$elemnr) + 0.5); #longitud media de los elementos genomicos
    if($assembly_dir){ #suponemos que nos dan la secuencia
      print $assembly_dir."\n";
      my $file = $assembly_dir.'/'.$chrom.'.fa';
      print "\n      ***                  Reading   Sequence                ***\n";
      ($seqst, $ID) = &GetFA($file); #seqst es toda la secuencia sin \n, ID es lo que hay detras de ">"

		  ($num_dinuc, $seqlen) = &GetSeqFeat($seqst);
      $seqlenbruto = length($seqst);
      $last_chunk = $seqlenbruto - $chrom_bed{$chrom}{'cEnd'}[$elemnr - 1]; #seqlenbruto?
      print "\nLastChunk:".$last_chunk."_\n";
      print "\nSeqlenbruto:".$seqlenbruto."_\n";
      print "\nUltimaCend:".$chrom_bed{$chrom}{'cEnd'}[$elemnr - 1]."_\n";
		}else{ #ESTIMACIONES
      #$seqlen = $chrom_bed{$chrom}{'cEnd'}[$elemnr - 1]; #estimacion = cEnd del ultimo elemento
      my $sum2 = $chrom_bed{$chrom}{'cStart'}[0];
      for (my $i = 1; $i < $elemnr; $i++){
        $sum2 += $chrom_bed{$chrom}{'cStart'}[$i] - $chrom_bed{$chrom}{'cEnd'}[$i-1];
      }
      $last_chunk = (int(($sum2/$elemnr) + 0.5)); # se hace un round()
      $seqlen = $chrom_bed{$chrom}{'cEnd'}[$elemnr - 1] + $last_chunk;
      $num_dinuc = $seqlen - 1;
      $seqlenbruto = $seqlen;

    }
    #push @{ $all_dist_n{$chrom} }, $last_chunk; #TODO: comprobar si tiene Ns antes de insertar
    my $nm = 0;
    if ($nBedFile){
      $nm = matchN($chrom, $chrom_bed{$chrom}{'cEnd'}[$elemnr - 1], $chrom_bed{$chrom}{'cEnd'}[$elemnr - 1] + $last_chunk);
    }
    #push @all_dist_n, $last_chunk if $nm <= $maxN;
    push @{ $all_dist_n{$chrom} }, $last_chunk if $nm <= $maxN;
    print "\nLastChunk:".$last_chunk."_\n";


		if($genome_intersec){ #genome
			$elemnr_genome += $elemnr;
			$length_genome += $num_dinuc;
      $elmean_genome += $elmean;
		}

		## Prob CpG ##
		#my $Ndach = $num_dinuc - $elemnr;
		#$prob = $elemnr/$Ndach;
    $prob = GetProb($num_dinuc, $elemnr, $elmean);

    #my @dist_n = [];
    #if($n_chrom_index{$chrom}[1] >= 0){
    #  @dist_n = @all_dist_n[$n_chrom_index{$chrom}[0] .. $n_chrom_index{$chrom}[1]];
    #}
    #push @dist_n, $last_chunk if $nm <= $maxN;

		#** Begin: single chrom INTERSEC
		if($chrom_intersec){
      #foreach (@{$all_dist_n{$chrom}}){
      #  print;
      #  print "\n";
      #}
      #die "...."
      print $#{$all_dist_n{$chrom}};

			print "      ***          Calculating  Chrom  Intersection          ***\n";
			my ($max, $min) = getMinMaxDistance($all_dist_n{$chrom}, $prob);
			$d = $max;
			$getd = "Chromosome Intersection";
			### get protoislands
			print "      ***      Detecting CpG clusters  (Chrom Intersec)      ***\n";
			my @protoislas = &GetProtoIslas($chrom_bed{$chrom}{'cStart'},$chrom_bed{$chrom}{'cEnd'});
			print "      ***       Calculating P-values  (Chrom Intersec)       ***\n";
			@protoislas = &CalcPvalNB(\@protoislas,$prob,$elmean);
			## Get Features like the obs/esp, clustering etc....
			#@protoislas = &GetCGI_features(\@protoislas);

			## Writing output
			$output .= $chrom."_chromIntersec_CpGcluster.txt";
      $elemmn = $elmean;
			&OUT_f(\@protoislas, $chrom);
			$output = $output_aux;
		}

		#** End: single chrom INTERSEC




		#** Begin: single chrom PERCENTILE
		if(@getd_hash > 0){
      foreach (@{$all_dist_n{$chrom}}){
        print;
        print "\n";
      }
			@dd = sort {$a <=> $b} @{ $all_dist_n{$chrom} }; #TODO: quitar variables globales como @dd

			foreach(@getd_hash){
				$getd = $_;
				print "      ***               Calculating Percentile $getd            ***\n";
				$d = &GetPerc(\@dd,$_);

				### get protoislands
				print "      ***             Detecting CpG clusters (p$getd)           ***\n";
				my @protoislas = &GetProtoIslas($chrom_bed{$chrom}{'cStart'},$chrom_bed{$chrom}{'cEnd'});



#open (P,">pinchado222.txt") or die "Can't open pinchado222.txt";
#foreach my $elem (@protoislas){
#  print P @{$elem};
#  print P "\n";
#  }
#  close (P);





				print "      ***              Calculating P-values  (p$getd)           ***\n";
				@protoislas = &CalcPvalNB(\@protoislas,$prob,$elmean);

        print "getting last features...\n";
				## Get Features like the obs/esp, clustering etc....
				#@protoislas = &GetCGI_features(\@protoislas);
        print "writing\n";
				## Writing output
				$output .= $chrom."_p".$getd."_CpGcluster.txt";
        $elemmn = $elmean;
				&OUT_f(\@protoislas, $chrom);
				$output = $output_aux;

			}
		}

		#** End: single chrom PERCENTILE

		push @dist_all,@{ $all_dist_n{$chrom} } if($genome_intersec); #genome
    #undef (@{ $all_dist_n{$chrom} }); # Podría dar problemas

	}


  #TODO: Ya no tiene sentido sacar la IG del bucle principal. Hay que ponerlo dentro para optimizar recursos.
  #Cambiar la estructura de $all_dist_n

	#** Begin: genome
	if($genome_intersec){
		print "\n\n      --- Genome:\n\n";
    #$elemnr = 0;
    #for my $chrom_key (keys %chrom_bed){
    #  $elemnr += $#{ $chrom_bed{$chrom_key}{'cStart'} } + 1;
    #}

		#my $Ndach_genome = $length_genome - $elemnr_genome;
		#my $prob_genome = $elemnr_genome/$Ndach_genome;
    $elmean_genome /= (scalar keys %chrom_bed);
    $prob = GetProb($length_genome, $elemnr_genome, $elmean_genome);
		my ($max, $min) = getMinMaxDistance(\@dist_all, $prob);

    $d = $max;
    $getd = "Genome Intersection";
		foreach my $chrom (keys %chrom_bed){
			print "\n      ".$chrom."\n";

			### get protoislands
			print "      ***      Detecting CpG clusters (Genome Intersec)      ***\n";
			#($seqst, $ID) = &GetFA($assembly_dir."/".$chrom.".fa");
			my @protoislas = &GetProtoIslas($chrom_bed{$chrom}{'cStart'},$chrom_bed{$chrom}{'cEnd'});
			print "      ***       Calculating P-values (Genome Intersec)       ***\n\n\n";
			@protoislas = &CalcPvalNB(\@protoislas,$prob,$elmean_genome);
			## Get Features like the obs/esp, clustering etc....
			#@protoislas = &GetCGI_features(\@protoislas);

			## Writing output
			$output .= $chrom."_genomeIntersec_CpGcluster.txt";
      $elemmn = $elmean_genome;
      $elemnr = $elemnr_genome;
			&OUT_f(\@protoislas, $chrom);
			$output = $output_aux;
		}
	}

#	closedir($DIR);
#}else{
#	die "Cannot open $assembly_dir\n"
#}

#** End: genome



####################################################################
#######   SUBFUNCTIONS   ###########################################
###################################################################

sub GetDefault{
  print "\n";
  print "---------------------------------------------------------------------------\n";
  print "---------------------------------------------------------------------------\n";
  print "---------                                                         ---------\n";
  print "---------   Laboratorio de Genomica Evolutiva y BioInformatica    ---------\n";
  print "---------    Universidad de Granada, Departamento de Genetica     ---------\n";
  print "---------                                                         ---------\n";
  print "---------               Web: http://bioinfo2.ugr.es               ---------\n";
  print "---------        CGI: http://bioinfo2.ugr.es/CpGcluster           ---------\n";
  print "---------                                                         ---------\n";
  print "---------              CpGcluster (2.0) 11/30/11                  ---------\n";
  print "---------                                                         ---------\n";
  print "---------------------------------------------------------------------------\n";
  print "---------------------------------------------------------------------------\n";
  print "\n";


	if($#ARGV < 2){
	    print "Example for the usage of CpGcluster:\n\n";
	    #print "perl CpGcluster.pl <assembly>  <d>  <P-value>\n\n";
      print "perl Program.pl <BED> <d> <P-value> [<assembly> [<N_BED> [<maxN>]]]";

	    print "\nassembly:   Directory containing sequence files in FASTA format\n";

	    print "\nd:          The threshold distance on basis of a given percentile.\n";
	    print "            For example: d=25 calculates the percentile 25 of the genomic\n";
	    print "            CpG distance distribution and takes this value as the threshold\n";
	    print "            distance\n";
	    print "            The recommended value is 50 (median distance)\n";
		print "            You can add multiple comma-separated percentile values, \"ci\"\n";
		print "            (chromosome intersection) or \"gi\" (genome intersection)\n";
		print "            Example: gi,25,60,ci,50\n";
	    print "\nP-value:    The maximal P-value under which a CpG cluster is considered as a\n";
	    print "            CpG island\n";
	    print "            The recommended limit is 1E-5\n\n";

	    die "\n";
	}

  #DEFAULT
  $assembly_dir = 0;
  $maxN = 0;
  $nBedFile = 0;


  my $i= 0;
  if(-e $ARGV[$i]){
    $bed_file = $ARGV[$i];
  }
  else{
    die "Cannot find the input file: $ARGV[$i]\n";
  }

  $i++;
    #ARGV1 - Percentile/ci/gi
	$getd = $ARGV[$i];
	foreach(split(',',$ARGV[$i])){
		if(/\D/){
			if(lc eq 'ci'){
				$chrom_intersec = 1;
			}elsif(lc eq 'gi'){
				$genome_intersec = 1;
			}
		}elsif($_ < 0 or $_ > 100){
			die "The Percentile must be between 0 and 100\n";
		}else{
			push @getd_hash,$_;
		}
	}

  $i++;
  	#ARGV2 - pLimit
    $plimit = $ARGV[$i];
	die "The maximal P-value you have choosen is higher than 1!\nPlease revise the order of the input parameters\n" if($plimit > 1);

  if($#ARGV > $i){
    $i++;
    #ARGV3 - assembly/directory
    if(-d $ARGV[$i]){
      $assembly_dir = $ARGV[$i];

    }else{
      die "Cannot find the input directory: $ARGV[$i]\n";
    }
    if($#ARGV > $i){
      $i++;
      if(-e $ARGV[$i]){
        $nBedFile = $ARGV[$i];
        getNbed($nBedFile);
      }else{
        die "Cannot find the inputo file: $ARGV[$i]";
      }

      if($#ARGV > $i){
        $i++;
        $maxN = $ARGV[$i];
        die "Incorrect value for maxN" if ($maxN < 0);

      }
    }
  }


  mkdir "./result", 0755;
	#my @f = split(/\//,$assembly_dir);
	#my @f1 = split(/\./,$f[$#f]);
	#$f[$#f] = "$f1[0]";
	#$output = join('/',@f);
	$output = "./result/";


}

# Read Sequence
sub GetFA{

  my $seqst_temp = "";
  open (I,$_[0]) or die "Can't open $_[0]";
  my $z = <I>;
  my $tes = substr($z,0,1);
  if($tes ne ">"){
    die "Sequence seems not to be in fasta format !!!!";
  }
  my @z = split(/\s+/," $z");
  $z = $z[1];
  $z=~ s/>//;
  $z =~ s/[\n\t\f\r\s]//g;
  my $ID_temp = $z;
  while($z = <I>){
    $z =~ s/[\n\t\f\r_0-9\s]//g;
    $seqst_temp .= $z;
  }
  close(I);

  return ($seqst_temp,$ID_temp);
}

sub GetBED{

  my $last_coord = 0;
  #my $i = 0;
  open (B,$bed_file) or die "Can't open $bed_file";
  #open (AUX, '>distN.txt');

  while ( my $line = <B> ) {
     my $rec = {};
     my @recsplit = split("\t", $line);
     my $chrom = $recsplit[0];

     #if(!exists($chrom_bed{$chrom})){ #Reinicializar para cuando cambie de cromosoma
     #  $last_coord = 0;
     #  $n_chrom_index{$chrom} = [$i, $i-1];
     #}

     push @{ $chrom_bed{$chrom}{'cStart'} }, $recsplit[1];
     push @{ $chrom_bed{$chrom}{'cEnd'} }, $recsplit[2];
     #push @{ $chrom_bed{$chrom}{'id'} }, $recsplit[3];
     #push @{ $chrom_bed{$chrom}{'score'} }, $recsplit[4];
     #push @{ $chrom_bed{$chrom}{'strand'} }, $recsplit[5];


     #TODO: comprobar si el trozo en cuestion contiene Ns
     my $nm = 0;
     if ($nBedFile){
      $nm = matchN($chrom, $last_coord, $recsplit[1]);
     }

     push @{ $all_dist_n{$chrom} }, ($recsplit[1] - $last_coord) if($nm <= $maxN);
     #print AUX ($recsplit[1] - $last_coord)."\n" if($nm <= $maxN);
     $last_coord = $recsplit[2];

     


     #if($nm <= $maxN){
     #  push @all_dist_n, ($recsplit[1] - $last_coord);
     #  $n_chrom_index{$chrom}[1] = $i; #indices [ ]
     #  $i++;
     #}
     #$last_coord = $recsplit[2];

  }
  #close(AUX);
}

sub OUT_f { #TODO: crear nuevo output
  my $c=0;

  open (OO,">$output") or die "Cannot not open $output";
  print OO "Chrom\tFrom\tTo\tLength\tCount\tPatDen\tPValue\tlogPValue\n";


  while($_[0]->[$c]){
	my $log_pvalue = ($_[0]->[$c]->[4] == 0 ? 0 : (log($_[0]->[$c]->[4])/log(10)));
	my $patden = (($_[0]->[$c]->[3])*$elemmn/$_[0]->[$c]->[2]);
	printf OO "%s\t%i\t%i\t%i\t%i\t%.3f\t%.2e\t%.2f\n", $_[1], $_[0]->[$c]->[0], $_[0]->[$c]->[1], $_[0]->[$c]->[2], $_[0]->[$c]->[3], $patden, $_[0]->[$c]->[4],$log_pvalue;
    $c++;
  }
  close(OO);
  open(O,">$output-log.txt") or die "can't open $output-log.txt";
  print O "Basic statistics of the input sequence: $_[1]"."_".$getd."\n";
  printf O "Length: %d\n",$seqlen;
  printf O "Length including Ns: %d\n",$seqlenbruto;
  #my $fg = $seqst =~ s/g/g/ig; #TODO: arreglar estas variables para Genome Intersec
  #my $fc = $seqst =~ s/c/c/ig;
  #my $fa = $seqst =~ s/a/a/ig;
  #my $ft = $seqst =~ s/t/t/ig;
  #printf O "GC content: %0.3f\n",100*($fg+$fc)/$seqlen;
  printf O "Number of elements in sequence: %d\n",$elemnr;
  printf O "Element mean length in sequence: %d\n",$elemmn;
  printf O "Probability to find a genomic element: %.4f\n\n",$prob;
  printf O "Calculated distance threshold: %i\n",$d;
  print O "Parameters used:\n";

  printf O "p-value threshold: $plimit\n";
  if($getd){
    print O "Distance threshold method: ";
    $_ = $getd;
    print O "percentile " if(/\d/);
    print O "$getd\n";
  }

  close(O);
}

sub OUT_f2 { #TODO: crear nuevo output
  my $c=0;

  open (OO,">$output") or die "could not open $output";
  print OO "CGI\tFrom\tTo\tLength\tCount\tOEratio\t%G+C\tPatDen\tPValue\tlogPValue\n";


  while($_[0]->[$c]){
  my $log_pvalue = ($_[0]->[$c]->[8] == 0 ? 0 : (log($_[0]->[$c]->[8])/log(10)));
  my $patden = ($_[0]->[$c]->[3]/$_[0]->[$c]->[2]);
  my $gc = ($_[0]->[$c]->[7]/$_[0]->[$c]->[2]);
  printf OO "%i\t%i\t%i\t%i\t%i\t%.3f\t%.2f\t%.3f\t%.2e\t%.2f\n",$c+1, $_[0]->[$c]->[0], $_[0]->[$c]->[1], $_[0]->[$c]->[2], $_[0]->[$c]->[3], $_[0]->[$c]->[4], $gc*100, $patden, $_[0]->[$c]->[8],$log_pvalue;
    $c++;
  }
  close(OO);
  open(O,">$output-log.txt") or die "can't open $output-log.txt";
  print O "Basic statistics of the input sequence: $ID\n";
  printf O "Length: %d\n",$seqlen;
  printf O "Length including Ns: %d\n",$seqlenbruto;
  my $fg = $seqst =~ s/g/g/ig; #TODO: arreglar estas variables para Genome Intersec
  my $fc = $seqst =~ s/c/c/ig;
  my $fa = $seqst =~ s/a/a/ig;
  my $ft = $seqst =~ s/t/t/ig;
  printf O "GC content: %0.3f\n",100*($fg+$fc)/$seqlen;
  printf O "Number of elements in sequence: %d\n",$elemnr;
  printf O "Probability to find a CpG: %.4f\n\n",$prob;
  print O "Parameters used:\n";

  printf O "p-value threshold: $plimit\n";
  if($getd){
    print O "Distance threshold method: ";
    $_ = $getd;
    print O "percentile " if(/\d/);
    print O "$getd\n";
  }

  close(O);
}
## Get CpG cluster
sub GetProtoIslas{

  my @cStart = @{$_[0]};
  my @cEnd = @{$_[1]};
  my @t;
  my ($start, $end, $elementnr);
  my $des = "no";
  print $#cStart."\n^^^^\n";
  for(my $i = 0; $i <= $#cStart - 1; $i++){

    my $dist = $cStart[$i+1]  - $cEnd[$i]; ## revisar coordenadas

    if($dist <= $d){
      if($des eq "no"){
	      $start = $cStart[$i];
        $elementnr = 1;
      }
      $end = $cEnd[$i+1]; ##  - 1?
      $elementnr++;
      $des = "yes";
    }
    elsif($dist > $d and $des eq "yes"){
      $des = "no";
      my @f = ($start, $end, $elementnr);
      push @t,\@f;
    }
  }
  if($des eq "yes"){
    my @f = ($start, $end, $elementnr);
    push @t,\@f;
  }
  return @t;
}

sub GetCGI_features{

  my @temp;
  my $c=0;
  while(defined($_[0]->[$c])){
    if($_[0]->[$c]->[4] < $plimit){
      #TODO: arreglar valores cpg
      my $len = $_[0]->[$c]->[1] - $_[0]->[$c]->[0] +1;
      (my $oe, my $cpgseq, my $gccont)= &CalcObsEsp($_[0]->[$c]->[0] -1,$len);
      my $coord1 = &GetCoord($cpgseq);
      (my $clust, my $meandist) = &GetClust($coord1);

      my $pval = $_[0]->[$c]->[4];
      $_[0]->[$c]->[4] = $oe;
      $_[0]->[$c]->[5] = $meandist;
      $_[0]->[$c]->[6] = $clust;
      $_[0]->[$c]->[7] = $gccont;
      $_[0]->[$c]->[8] = $pval;
      my @t = @{$_[0]->[$c]};
      push @temp,\@t;
    }
    $c++;
  }
  return @temp;
}

sub CalcObsEsp{

  my $cpgseq = substr ($seqst,$_[0],$_[1]);
  my $fc = $cpgseq =~ s/c/c/ig;
  my $fg = $cpgseq =~ s/g/g/ig;
  my $CGICpG = $cpgseq =~ s/CG/CG/ig;
  my $e = $fc*$fg;
  my $oet = $CGICpG*length($cpgseq)/$e;
  return ($oet,$cpgseq,$fc+$fg);
}

sub GetSeqFeat{

  my $n = $_[0];
  #my $elemnr = $n =~ s/CG/Cg/ig;
  my $ln = length($n);
  my $NN = $n =~ s/N/N/ig;
  my $seqlen = $ln-$NN;
  my $num_dinuc = 0;


  my $elem;
  my $limit = $ln - 1;
  for(my $i = 0; $i < $limit; $i++){
  	$elem = substr($n,$i,2);
  	if($elem !~ m/N/i){
  		$num_dinuc++;
  	}
  }
  return ($num_dinuc,$seqlen);
}

sub GetProb{
  my ($length, $elemnr, $elmean) = @_;
  return (($elemnr * $elmean)/$length);
}

sub GetCoords{

  my $n = $_[0];

  $n.="j";
  my @f =split(/$_[1]/i,$n);
  my @t;

  my $lencount = 0;
  for(my $i = 0; $i < $#f; $i++){
    $lencount += length($f[$i]);
    $t[$i] = $lencount + 1;
    $lencount+=2;
    my $nnr = $f[$i] =~ s/n/n/ig;
    if($nnr <= $maxN){
      push @dist_n,length($f[$i])+1;
    }
  }
  return @t;
}

sub GetPerc{

  my @t = @{$_[0]};
  my $totnr = @t;
  for(my $i = 0; $i <= $#t; $i++){
    if(100*$i/$totnr >= $_[1]){
      return $t[$i];
    }
  }
  return $t[$#t];
}

sub getNbed{
  open(N, '<'.$_[0]) or die "Cannot open file input: $_[0]\n";
    while ( my $line = <N> ) {
      my $rec = {};
      my @recsplit = split("\t", $line);
      my $chrom = $recsplit[0];

      if(!exists($n_coords{$chrom})){
        $n_coords{$chrom} = [];
      }

      push $n_coords{$chrom}, [$recsplit[1], $recsplit[2]];
  }

  close(N);

}

sub matchN{
  my ($chrom, $cStart, $cEnd) = @_;

  foreach my $elem (@{$n_coords{$chrom}}){


    if(($cStart >= $elem->[0] && $cStart <= $elem->[1]) || ($cEnd >= $elem->[0] && $cEnd <= $elem->[1])){
print "o\n";
      return ($elem->[1] - $elem->[0]); #hay N
    }elsif($elem->[0] > $cEnd){
      last;

    }
  }
  print "q";
  return 0; #No hay N
}


##########################################################################
############## SUBFUNCTIONS for P-value calculations  ####################
##########################################################################

######################################################################
## *** Calculates the negative binomial distribution

sub CalcPvalNB{

  my ($pval, @temp);
  my $c=0;
  my @islas = @{$_[0]};
  my $prob_ge = $_[1];
  my $elmean = $_[2];
print "debug1\n";

  for(my $i = 0; $i <= $#islas; $i++){
print $i;
    my $l = $islas[$i]->[1] - $islas[$i]->[0];
          #my $str = substr ($seqst,$islas[$i]->[0]-1,$l);
          #my $cpg = $str =~ s/cg/cg/ig;
    my $gel = $islas[$i]->[2];
          #my $pval = &GetNB({{{$l-(2*$cpg}}}),$cpg-1,$_[1]);
    my $nf = int(($l-($elmean*$gel))/$elmean) + 1; #ceiling
    my $pval = &GetNB($nf,$gel-1,$prob_ge); #Mirar que el numero de no-elementos sea positivo
    $pval = sprintf("%.5e",$pval);
    my @t = ($islas[$i]->[0],$islas[$i]->[1],$l,$gel,$pval);
    push @temp, \@t;
    $c++;
  }

  print "\ndebug2\n";
  return @temp;
}
sub GetNB{
  my $pval = 0;
  for(my $j = 0; $j <= $_[0]; $j++){
    my $ptemp = &FactorialNB($j,$_[1]) + $_[1]*log($_[2]) + $j*log(1.0-$_[2]);
    $ptemp = exp($ptemp);
    $pval += $ptemp;
    }
  return $pval;
}


sub FactorialNB{
  my $stop = $_[0]+$_[1]-1;
  my $l1 = 0;
  my $l2 = 0;
  for(my $i = $_[0]+1;$i <= $stop; $i++){
    $l1 +=log($i);
  }
  for(my $i = 1;$i <= $_[1]-1; $i++){
    $l2 +=log($i);
  }
  return $l1-$l2;
}
##########################################################################
############## SUBFUNCTIONS for Clustering Calculation  ####################
##########################################################################

sub GetClust{

  my @t = @{$_[0]};
  my $dist = &GetDist($_[0]);
  my $mean = &Normalize($dist);
  return (1, $mean);
}

sub Normalize{

  my @d = @{$_[0]};
  my $tot;
  for(my $i = 0; $i <= $#d; $i++){
    $tot+=$d[$i];
  }
  my $mean = $tot/@d;
  return $mean;
}

sub GetDist{

  my @dist;
  my @d = @{$_[0]};
  for(my $i = 0; $i < $#d; $i++){
    my @f = split (/\-/,$d[$i]);
    my @s = split (/\-/,$d[$i+1]);
    my $dist = $s[0]-$f[1];
    push @dist,$dist;
  }
  return \@dist;
}

sub GetCoord{

  my @coord;
  my @c = split (//,$_[0]);
  for(my $i = 0; $i < $#c; $i++){
    if($c[$i] =~ /[cC]/ and $c[$i+1] =~ /[gG]/){
      my $str = $i.'-'.eval($i+1);
      push @coord,$str;
    }
  }
  return \@coord;
}


sub getMinMaxDistance{
	my @distances = @{$_[0]};
	my $prob = $_[1];

	my @distCount = ();
	my $maxDist = 0;
	my $nrDistances = 0;
	my $stop = @distances;

	for(my $i = 0; $i < $stop; $i++){
		$nrDistances++;
		my $dist = $distances[$i];

		if(defined $distCount[$dist]){
			$distCount[$dist]++;
		}else{
			$distCount[$dist] = 1;
		}

		if($dist > $maxDist){
			$maxDist = $dist;
		}
	}

	my $obsCum = 0;
	my $teoCum = 0;
	my @dif = ();
print "nrDistances:".$nrDistances."\n";
	for(my $i = 1; $i <= $maxDist; $i++){
		if(defined $distCount[$i]){
			$obsCum += $distCount[$i]/$nrDistances;
		}
		$teoCum += $prob * (1 - $prob)**($i-1);
		$dif[$i] = ($obsCum - $teoCum);
	}

	my $max = -1;
	my $min = 1;
	my $maxD = 0;
	my $minD = 0;
	for(my $i = 1; $i <= $maxDist; $i++){
		if(defined $dif[$i]){
			my $difT = $dif[$i];
			if($difT > $max){
				$max = $difT;
				$maxD = $i;
			}
			elsif($difT < $min){
				$min = $difT;
				$minD = $i;
			}
		}
	}
  print $maxD."   ---   ".$minD."\n";
	return ($maxD, $minD);
}
#! /usr/bin perl

# Codes from mastering algorithms with perl
# choose($n, $k) is the number of ways to choose $k elements from a set
# of $n elements, when the order of selection is irrelevant.
#

sub choose {
  my ($n, $k) = @_;
  my ($result, $j) = (1, 1);
  return 0 if $k > $n || $k < 0;
  $k = ($n - $k) if ($n - $k) < $k;
  while ( $j <= $k ) {
    $result *= $n--;
    $result /= $j++;
  }
  return $result;
}


sub hypergeometric {
  
  
  my ($x, $k, $m, $n) = @_;
  return unless $m > 0 && $m == int($m) && $n > 0 && $n == int($n) &&
    $k > 0 && $k <= $m + $n;
  return 0 unless $x <= $k && $x == int($x);
  return choose($m, $x) * choose($n, $k - $x) / choose($m + $n, $k);
}
sub hypergeometric_expected { $_[0] * $_[1] / ($_[1] + $_[2]) }
sub hypergeometric_variance {
  my ($k, $m, $n) = @_;
  return $m * $n * $k * ($m + $n - $k) / (($m + $n) ** 2) /
    ($m + $n - 1);
}

sub sign_significance {
  my ($phased, $totalreads, $phased_pos, $nonphased_pos) = @_;
  my $confidence;
  foreach ($phased..$phased_pos) {
    $confidence += hypergeometric ($_, $totalreads, $phased_pos, $nonphased_pos);
  }
  return $confidence;
}




# The scan method for phased loci (Pvalue) and phasing score is different.
#
# For phased loci, each mapped read will initialize a 9*21bp scan window, 
# and shift with 3*21bp, and extend if the region meet the criteria (>=10 
# unique total reads, >=3 phased reads, >=50% total reads are 21bp). 
#
# For phasing score, each cordinate of chromosome will initialize a scan
# window of fixed length of 9*21bp.

# Input mapping file should follow format below (separated by tabular):
#
# Chr, start, end, readid, strand, sequence
#
# readid example: seq123-50, 50 is read abundance. 
print "Start at ".localtime."\n";
open(MAP, "<", $ARGV[0]) || die "Cannot open mapping file $ARGV[0]: $!";
open(OUT, ">", $ARGV[1]) || die "Cannot open $ARGV[1]: $!";
my $RegisterLen = $ARGV[2];
my %readsinfor;
my @loci;
my @phased;


while(<MAP>){
  chomp;
  my ($chr, $start, $end, $seqid, $strand, $sequence) = split(/\t/, $_);
  $readsinfor{$chr}{$strand}{$start}{$end} = $_;
}
close MAP;

# Start window scan for phased loci discovery
foreach my $key1 (keys %readsinfor){
  foreach my $key2 (keys %{$readsinfor{$key1}}){
    foreach my $key3 (keys %{$readsinfor{$key1}{$key2}}){
      foreach my $key4 (keys %{$readsinfor{$key1}{$key2}{$key3}}){
	my ($chr, $start, $end, $seqid, $strand, $sequence) = split(/\t/, $readsinfor{$key1}{$key2}{$key3}{$key4});
	# print "$chr,$start, $end, $seqid, $strand, $sequence\n";
	my ($loci_chr, $loci_start, $loci_end, $loci_phased, $loci_total, $loci_nt21);
	my ($win_start, $win_end, $totalreads, $phased, $nt21_reads, $extend, %phasedreads);
	
	if($end-$start+1== $RegisterLen){
	  $loci_chr = $chr;
	  $loci_start = $start;
	  $loci_end = $loci_start + 9*$RegisterLen -1;
	  $extend = 0;
	  %phasedreads=();
	  
	  do{    
	    $phased = 0;
	    $totalreads = 0;
	    $nt21_reads = 0;
	    
	    # define window of sense strand
	    if($strand eq "+"){
	      $win_start = $start + $extend;
	      $win_end = $win_start + 9*$RegisterLen-1 + $extend;
	    }else{
	      $win_start = $start + 2 + $extend;
	      $win_end = $win_start + 9*$RegisterLen-1 + $extend;
	    }
	    # calculate total reads on sense strand
	    for(my $i=$win_start; $i<=$win_end; $i++){
	      if(exists $readsinfor{$chr}{"+"}{$i}){
		my $size = keys %{$readsinfor{$chr}{"+"}{$i}};
		$totalreads = $totalreads + $size;
		foreach my $n (keys %{$readsinfor{$chr}{"+"}{$i}}){
		  if($n - $i + 1 == $RegisterLen){
		    $nt21_reads++;
		  }
		}
	      }
	    }
	# calculate phased reads on sense strand
	    for(my $i=$win_start; $i<=$win_end; $i=$i+$RegisterLen){
	      if(exists $readsinfor{$chr}{"+"}{$i}){
		foreach my $n (keys %{$readsinfor{$chr}{"+"}{$i}}){
		  if($n - $i + 1 == $RegisterLen){
		    $phased++;
		  }
		}
	      }
	    }
	    # define window of antisense strand
	    if($strand eq "+"){	    
	      $win_start = $start - 2 + $extend;
	      $win_end = $win_start + 9*$RegisterLen-1 + $extend;
	    }else{
	      $win_start = $start + $extend;
	      $win_end = $win_start + 9*$RegisterLen-1 + $extend;
	    }
	    # calculate total reads on antisense strand
	    for(my $i=$win_start; $i<=$win_end; $i++){
	      if(exists $readsinfor{$chr}{"-"}{$i}){
		my $size = keys %{$readsinfor{$chr}{"-"}{$i}};
		$totalreads = $totalreads + $size;
		foreach my $n (keys %{$readsinfor{$chr}{"-"}{$i}}){
		  if($n - $i + 1 == $RegisterLen){
		    $nt21_reads++;
		  }
		}
	      }
	    }
	    # calculate phased reads on antisense strand
	    for(my $i=$win_start; $i<=$win_end; $i=$i+$RegisterLen){
	      if(exists $readsinfor{$chr}{"-"}{$i}){
		foreach my $n (keys %{$readsinfor{$chr}{"-"}{$i}}){
		  if($n - $i + 1 == $RegisterLen){
		    $phased++;
		  }
		}
	      }
	    }
	    $loci_end += $extend;
	    $extend += 3*$RegisterLen;
	  }while($phased >= 3 and $totalreads >=10 and $nt21_reads/$totalreads >=0.5);
	  
	  # calculate for loci (extended window)
	  if($strand eq "+"){
	    $win_start = $loci_start;
	    $win_end = $loci_end;
	  }else{
	    $win_start = $loci_start + 2;
	    $win_end = $loci_end +2;
	  }
	  for(my $i=$win_start; $i<=$win_end; $i++){
	    if(exists $readsinfor{$chr}{"+"}{$i}){
	      my $size = keys %{$readsinfor{$chr}{"+"}{$i}};
	      $loci_total += $size;
	    }
	  }
	  for(my $i=$win_start; $i<=$win_end; $i=$i+$RegisterLen){
	    if(exists $readsinfor{$chr}{"+"}{$i}){
	      foreach my $n (keys %{$readsinfor{$chr}{"+"}{$i}}){
		if($n - $i + 1 == $RegisterLen){
		  $loci_phased ++;
		  $phasedreads{$readsinfor{$chr}{"+"}{$i}{$n}} = $reads{$readsinfor{$chr}{"+"}{$i}{$n}};
		}
	      }
	    }
	  }
	  if($strand eq "+"){
	    $win_start = $loci_start - 2;
	    $win_end = $loci_end - 2;
	  }else{
	    $win_start = $loci_start;
	    $win_end = $loci_end;
	  }
	  for(my $i=$win_start; $i<=$win_end; $i++){
	    if(exists $readsinfor{$chr}{"-"}{$i}){
	      my $size = keys %{$readsinfor{$chr}{"-"}{$i}};
	      $loci_total += $size;
	    }
	  }
	  for(my $i=$win_start; $i<=$win_end; $i=$i+$RegisterLen){
	    if(exists $readsinfor{$chr}{"-"}{$i}){
	      foreach my $n (keys %{$readsinfor{$chr}{"-"}{$i}}){
		if($n - $i + 1 == $RegisterLen){
		  $loci_phased ++;
		  $phasedreads{$readsinfor{$chr}{"-"}{$i}{$n}}=$reads{$readsinfor{$chr}{"+"}{$i}{$n}};
		}
	      }
	    }
	  }

	  if($loci_phased >=3 and $loci_total >=10){
	    my $pvalue = sign_significance($loci_phased, $loci_total, ($loci_end-$loci_start+1)/$RegisterLen, ($RegisterLen-1)*($loci_end-$loci_start+1)/$RegisterLen);
	    if($strand eq "+"){
	      push @loci, [ $loci_chr, $loci_start, $loci_end, $loci_phased, $loci_total, $pvalue ];
	    }else{
	      push @loci, [ $loci_chr, $loci_start+2, $loci_end+2, $loci_phased, $loci_total, $pvalue ];
	    }
	  }
	}
      }
    }
  }
}

# Sort locis
my @sorted_loci = sort {
    $a->[0] cmp $b->[0] 
	or
	$a->[1] <=> $b->[1]
	or
	$a->[2] <=> $b->[2]
} @loci;


#use Data::Dumper qw(Dumper);
#print Dumper \@sorted_loci;

# Merge locis
my @merged_loci;


for(my $i=0; $i<=$#sorted_loci; $i++){
  my ($chr, $start, $end) = ($sorted_loci[$i][0], $sorted_loci[$i][1], $sorted_loci[$i][2]);
  if($i==0){
    push @merged_loci, [$chr, $start, $end];
  }else{
    my $switch=1; # push switch
    for(my $z=0; $z<=$#merged_loci; $z++){
      if($chr eq $merged_loci[$z][0] and $start <= $merged_loci[$z][2] and $end > $merged_loci[$z][2] and ($start-$merged_loci[$z][1])%$RegisterLen==0){
	# make sure merged locis are in the same phasing increment
	$merged_loci[$z][2] = $end;
	$switch=0;
      }elsif($chr eq $merged_loci[$z][0] and $start < $merged_loci[$z][2] and $end <= $merged_loci[$z][2] and ($start-$merged_loci[$z][1])%$RegisterLen==0){
	$switch=0;
      }else{
	# do nothing;
      }
    }
    if($switch==1){
      push @merged_loci, [$chr, $start, $end];
    }
  }
}

# calculate for merged loci
my @lociphased;
    
for(my $i=0; $i<=$#merged_loci; $i++){
    my ($chr, $start, $end) = ($merged_loci[$i][0], $merged_loci[$i][1], $merged_loci[$i][2]);
    # Only 21nt reads are involved into calculation of pvalue, but total reads will also be recorded
    my ($phasedreads, $total, $total21nt);
    my ($phasedreads_abundance, $total_abundance, $total21nt_abundance);
    $lociphased[$i] = [];
    for(my $k=$start; $k<=$end; $k++){
      if(exists $readsinfor{$chr}{"+"}{$k}){
	my $size = keys %{$readsinfor{$chr}{"+"}{$k}};
	$total += $size;
	foreach my $n (keys %{$readsinfor{$chr}{"+"}{$k}}){
	  my @data = split(/\t/, $readsinfor{$chr}{"+"}{$k}{$n});
	  my ($seqname, $count) = split(/-/, $data[3]);
	  $total_abundance += $count;
	  if ($n - $k +1 == $RegisterLen){
	    $total21nt ++;
	    $total21nt_abundance += $count;
	  }
	}
      }
    }
    for(my $k=$start; $k<=$end; $k=$k+$RegisterLen){
      if(exists $readsinfor{$chr}{"+"}{$k}){
	foreach my $n (keys %{$readsinfor{$chr}{"+"}{$k}}){
	  my @data = split(/\t/, $readsinfor{$chr}{"+"}{$k}{$n});
	  my ($seqname, $count) = split(/-/, $data[3]);
	  if($n - $k + 1 == $RegisterLen){
	    $phasedreads ++;
	    $phasedreads_abundance += $count;
	    push @{$phased[$i]}, $readsinfor{$chr}{"+"}{$k}{$n};
	  }
	}
      }
    }
    for(my $k=$start-2; $k<=$end-2; $k++){
      if(exists $readsinfor{$chr}{"-"}{$k}){
	my $size = keys %{$readsinfor{$chr}{"-"}{$k}};
	$total += $size;
	foreach my $n (keys %{$readsinfor{$chr}{"-"}{$k}}){
	  my @data = split(/\t/, $readsinfor{$chr}{"-"}{$k}{$n});
	  my ($seqname, $count) = split(/-/, $data[3]);
	  $total_abundance += $count;
	  if ($n - $k +1 == $RegisterLen){
	    $total21nt ++;
	    $total21nt_abundance += $count;
	  }
	}
      }
    }
    for(my $k=$start-2; $k<=$end-2; $k=$k+$RegisterLen){
      if(exists $readsinfor{$chr}{"-"}{$k}){
	foreach my $n (keys %{$readsinfor{$chr}{"-"}{$k}}){
	  my @data = split(/\t/, $readsinfor{$chr}{"-"}{$k}{$n});
	  my ($seqname, $count) = split(/-/, $data[3]);
	  if($n - $k + 1 == $RegisterLen){
	    $phasedreads ++;
	    $phasedreads_abundance += $count;
	    push @{$phased[$i]}, $readsinfor{$chr}{"-"}{$k}{$n};
	  }
	}
      }
    }
    my $pvalue = sign_significance($phasedreads, $total21nt, ($end-$start+1)*2/$RegisterLen, ($RegisterLen-1)*($end-$start+1)*2/$RegisterLen);
    push @{$merged[$i]}, ($phasedreads, $total21nt, $total, $pvalue);
    print OUT "$chr\t$start\t$end\t$phasedreads\t$total21nt\t$total\t$phasedreads_abundance\t$total21nt_abundance\t$total_abundance\t$pvalue\n";
    print OUT "-----\n";
    for(my $j=0; $j<=@{$phased[$i]}; $j++){
	print OUT $phased[$i][$j]."\n";
    }
    print OUT "//\n";
}

close OUT;
print "End at ".localtime."\n";


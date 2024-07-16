#! /usr/bin/perl -w

use strict;

my ($line, @data);
my ($allele, $name, $query, $target, $Qrange1, $Qrange2, $Trange1, $Trange2);
my ($point, $sub, @data1, @data2, $i);
my (%sr, %pe);

my $sample   = $ARGV[0];
my $retro    = $ARGV[1];
my $path     = $ARGV[2];
my $genotype = $ARGV[3];
my $cord1    = $ARGV[4];
my $cord2    = $ARGV[5];
my $exo_file = $path.'/'.$retro.'/'.$sample.'.L1HS.alignout';
my $ou1_file = $path.'/'.$retro.'/'.$sample.'.discover';
my $ou2_file = $path.'/'.$retro.'/'.$sample.'.sr.discover';
my $sr_file  = $path.'/'.$retro.'/'.$sample.'.sr.L1HS.discover';
my $pe_file  = $path.'/'.$retro.'/'.$sample.'.L1HS.discover';

(open (EXO, "<$exo_file")) || die "cannot open the exo file\n";
(open (OU1, ">>$ou1_file")) || die "cannot open the out file\n";
(open (OU2, ">>$ou2_file")) || die "cannot open the out file\n";
(open (SR , "<$sr_file ")) || die "cannot open the sr file\n";
(open (PE , "<$pe_file ")) || die "cannot open the pe file\n";

while ($line = <SR>)
   {
    chomp($line);
    @data = split("\t", $line);
    $sr{'sr_'.$data[3].'_'.$data[4].'_'.$data[8]} = $line;
   }

while ($line =<PE>)
   {
    chomp($line);
    @data = split("\t", $line);
    $pe{'pe_'.$data[3].'_'.$data[4].'_'.$data[7]} = $line;
   }

while ($line = <EXO>)
   {
    chomp($line);
    if ($line =~ /Query\:/)
       {
        @data = split(" ", $line);
        $name = $data[1];
        <EXO>;
        <EXO>;
        <EXO>;
        chomp($line = <EXO>);
        @data = split(" ", $line);
        $Qrange1 = $data[2];
        $Qrange2 = $data[4];
        chomp($line = <EXO>);
        @data = split(" ", $line);
        $Trange1 = $data[2];
        $Trange2 = $data[4];
        $query = $target = '';

        <EXO>;
        while ( ($line = <EXO>) && ($line !~ /vulgar/))
            {
             chomp($line);
             @data = split(" ", $line);
             $query .= $data[2];
             <EXO>;
             chomp($line = <EXO>);
             @data = split(" ", $line);
             $target .= $data[2];
             <EXO>;
            }
       if ( ($Trange1 < $cord1 && $Trange2 >= $cord2) || ($Trange1 >= $cord2 && $Trange2 < $cord1) )
          {
           $allele = &parse_genotype($query, $target, $Trange1, $Trange2, $cord1, $cord2);
       
           # print "$Trange1\t$Trange2\t$allele\n";
           ### print the updated discover ###
           if ($allele eq $genotype) # || $allele eq "ACG")
             {
              if ($name =~ /pe/)
                 {
                  print OU1 "$pe{$name}\t$allele\n";
                 }
              elsif ($name =~ /sr/)
                 {
                  print OU2 "$sr{$name}\t$allele\n";
                 }
             }
          }
      }
   }
       
sub parse_genotype
   {
    my $qseq  = shift;
    my $tseq  = shift;
    my $tnum1 = shift;
    my $tnum2 = shift;
    my $c1    = shift;
    my $c2    = shift;
    my $allele;

    my @data1 = split("", $tseq);
    my @data2 = split("", $qseq);

    if ($tnum1 < $tnum2)
       {
        $point = $tnum1;
        $sub = '';
        for ($i=0; $i<=$#data1; $i++)
           {
            if ($data1[$i] eq '-')
               {
                $point--;
               }
            $point++;
            if ($point >= $c1 && $point <= $c2)
               {
                $sub .= $data2[$i];
               }
            last if ($point > $c2 + 5);
           }
         $allele = $sub; 
       }
    elsif ($tnum1 > $tnum2) 
       {
        $point = $tnum1;
        $sub = '';
        for ($i=0; $i<=$#data1; $i++)
           {
            if ($data1[$i] eq '-')
               {
                $point++;
               }
            $point--;
            if ($point >= ($c1-1) && $point <= ($c2-1))
               {
                $sub .= $data2[$i];
               } 
            last if ($point < $c1-5);
           }
        $allele = &rc($sub);
       }
    return ($allele);
   }
    
sub rc 
  { ### reverse complement ###
   my $seq = shift;
   my %RCbase = (
      'A' => 'T',
      'T' => 'A',
      'G' => 'C',
      'C' => 'G',
      '-' => '-',
      'N' => 'N',
     );
   my @bases = split("", $seq);
   my @output;
   my ($seqout, $i);

   for ($i=0; $i<=$#bases; $i++)
     {
      $output[$i] = $RCbase{$bases[length($seq)-$i-1]};
     }
   $seqout=join("", @output);
   #print "$seq\n$seqout\n";
   return($seqout);
  }


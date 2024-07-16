#! /usr/bin/perl -w

### filter anchor ends    ###
### filter mapping score  ###
### filter mapping length ### 
### filter on avg seq depth ###

use strict;

my $subject   = $ARGV[0];
my $path      = $ARGV[1];
my $retro     = $ARGV[2];
my $pe_file   = $path.'/'.$retro.'/'.$subject.'.discover';
my $ou_file   = $path.'/'.$retro.'/'.$subject.'.alignfilter.discover';
my $o2_file   = $path.'/'.$retro.'/'.$subject.'.transduct.signals';
my $fa_file   = $path.'/'.$retro.'/seg.pe.fasta';
my $dep1_file = $path.'/'.$retro.'/'.$subject.'.bed';
my $dep2_file = $path.'/'.$retro.'/'.$subject.'.depth';
my ($c1, $c2, $c3, $direction, $min, $max, $size, $tag);
my %size = (
   'AluYa5' => 282,
   'AluYk13' => 282,
   'AluYb8' => 289,
   'AluYb9' => 289,
   'AluYc1' => 282,
   'AluYa5a2' => 282,
   'L1HS' => 6064,
   'HERVK' => 7536,
   'HERVH' => 7713,
   'SVA_E' => 1382,
   'SVA_F' => 1375,
   'L1' => 6064,
   'Alu' => 289,
   'SVA' => 1382,
    );

my ($CUT1, $CUT2, $CUT3, $CUT4);
my ($sim_ALU, $sim_L1, $sim_other);
my $FILTER = $ARGV[3];                ### 0=no_filter, 1=pre_filter, 2=all_filter ###
my $strand = ($ARGV[4] == 1) ? 1 : -1; ### 1: + strand ; 0: - strand ###
if ($FILTER == 2)
   {
    $sim_ALU = 92.67;
    $sim_L1  = 98;
    $sim_other = 95;
    $CUT1 = 80;
    $CUT2 = 40;
    $CUT3 = 1500; ### longest fragment in sequencing ###
    $CUT4 = 200;  ### size cutoff for too short L1 insertions ###
   }
else
   {
    $sim_ALU = 85;
    $sim_L1  = 85;
    $sim_other = 85;
    $CUT1 = 10;
    $CUT2 = 10;
    $CUT3 = 1500; ### longest fragment in sequencing ###
    $CUT4 = 200;  ### size cutoff for too short L1 insertions ###
   }

### Anchor ###
#### mapped > 85% ###
##### no alternative alignment ###
##### mismatch <= 2 in 50 bp ###
my $MM = 2; ## in read  of 50 bp ###

(open (PE,   "<$pe_file  ")) || die "cannot open the pe file\n";
(open (DEP1, "<$dep1_file")) || die "cannot open the dep1 file\n";
(open (DEP2, "<$dep2_file")) || die "cannot open the dep2 file\n";
(open (OU,   ">$ou_file ")) || die "cannot open the ou file\n";
(open (FA,   ">$fa_file ")) || die "cannot open the fa file\n";
(open (O2,   ">$o2_file ")) || die "cannot open the o2 file\n";

my ($line, @data);
my ($te, $i, $j, $k, $id, $len, $read, $r1, $seq, $seq1, $seq2, $pass);
my ($a_i, $a_size, $a_map, $a_MM, $med, $mad, $med_x, $mad_x, $med_y, $mad_y);
my (@temp1, %PE_cord, %depth, @dep_X, @dep_Y, @dep_all, $PE_depth);

### read the sequencing depth ###
while ($line = <DEP1>)
   {
    chomp($line);
    @data = split("\t", $line);
    @temp1 = split('__', $data[3]);

    if ($line =~ /PE/)
       {
        $PE_cord{$temp1[1]} = $data[0]."\t".$data[1]."\t".$data[2];
       }
   }

$i = 0;
$j = 0;
$k = 0;
while ($line = <DEP2>)
   {
    chomp($line);
    @data = split("\t", $line);
    $depth{$data[0]."\t".$data[1]."\t".$data[2]} = $data[3] / ($data[2] - $data[1]);
    if ($data[0] =~ /X/)
       {
        $dep_X[$i++] = $data[3] / ($data[2] - $data[1]);
       }
    elsif ($data[0] =~ /Y/)
       {
        $dep_Y[$j++] = $data[3] / ($data[2] - $data[1]);
       }
    else
       {
        $dep_all[$k++] = $data[3] / ($data[2] - $data[1]);
       }    
   }

$med = &median(@dep_all);
$mad = &mad(@dep_all);
$med_x = &median(@dep_X);
$mad_x = &mad(@dep_X);
$med_y = &median(@dep_Y);
$mad_y = &mad(@dep_Y);
print "$med\t$mad\t$med_x\t$mad_x\t$med_y\t$mad_y\n";

if ($FILTER == 2)
   {
    $sim_ALU = 92.67;
    $sim_L1  = 98;
    $sim_other = 95;
    $CUT1 = 80;
    $CUT2 = 40;
    $CUT3 = 1500; ### longest fragment in sequencing ###
    $CUT4 = 200;  ### size cutoff for too short L1 insertions ###
   }
else
   {
    $sim_ALU = 85;
    $sim_L1  = 85;
    $sim_other = 85;
    $CUT1 = 10;
    $CUT2 = 10;
    $CUT3 = 1500; ### longest fragment in sequencing ###
    $CUT4 = 200;  ### size cutoff for too short L1 insertions ###
   }
### Anchor ###
### mapped > 85% ###
#### no alternative alignment ###
#### mismatch <= 2 in 50 bp ###

while ($line = <PE>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $c1 = ($data[5] eq '+') ? 1 : -1;
    $c2 = ($data[6] eq '+') ? 1 : -1;
    $c3 = ($data[15] eq '+') ? 1 : -1;
    $direction = -1 * $c1 * $c2 * $c3;
    $tag = 'OK';
    next if ($direction ne $strand);
    $min = ($data[11] < $data[10]) ? $data[11] : $data[10];
    $max = ($data[11] > $data[10]) ? $data[11] : $data[10];

    if ($data[5] eq '+')
       { ### upstream read ###
        if ($direction > 0)
           { ### positive strand 5->3 ###
            if (($data[3] !~ /Alu/) && ($min > ($size{$data[3]}-$CUT4) ))
               {
                $tag = 'short';
               }            
           }
        else
           { ### negative strand 3->5 ###
            if ($max  < ($size{$data[3]} - $CUT3))
               {
                $tag = 'no3end';
               }
           }
       }
    else
       { ### downstream read ###
        if ($direction > 0)
           { ### positive strand 5->3 ###
            if ($max  < ($size{$data[3]} - $CUT3))
               {
                $tag = 'no3end';
               }
           }
        else
           { ### negative strand 3->5 ###
            if (($data[3] !~ /Alu/) && ($min > ($size{$data[3]}-$CUT4) ))
               {
                $tag = 'short';
               }
           }
       }

    next if ( ($tag eq 'no3end') && ($FILTER == 2));
    next if ( ($tag eq 'short') && ($FILTER == 2));

    $te   = $data[3];
    $id   = $data[4];
    $len  = abs($data[9] - $data[8]);  ### frag length mappable to TE ###
    $size = length($data[12]);      ### sequence length ###
 
    ### no XA ###
    next if ( ($line !~ /NULL/) && ($FILTER==2));    

    $a_size = $a_map = 0;
    while ($data[16] =~ /(\d+)/g)
       {
        $a_size += $1;                 ### anchor size, all parts ###
       }
    while ($data[16] =~ /(\d+)M/g)
       {
        $a_map += $1;                  ### anchor, all mappable size ###
       }

    ### enough mapped length on anchor reads ###
    next if ( (($a_map / $a_size) < (85/100)) && ($FILTER == 2));
    ### not too many mismatches ###
    next if ( ($data[18] > $MM * $a_size / 50) && ($FILTER==2) );
    
    ### remove TE alignment with too low similarities ##
    $read = $id.'###'.$te.'###'.$data[7].'###'.$data[8].'###'.$data[9];
    $seq  = substr($data[12], $data[8], $len);
    print FA ">$read\n$seq\n";

    if ($te =~ "Alu")
       {
        next if ($data[7] <= $sim_ALU);
       }
    elsif ($te =~ "L1")
       {
        next if ($data[7] <= $sim_L1);
       }
    else
       {
        next if ($data[7] <= $sim_other);
       }
   next if ((!exists($PE_cord{$r1.':'.$data[4]})) || (!exists($depth{$PE_cord{$r1.':'.$data[4]}})) );

   if ($data[0] =~ /X/)
      {
       $PE_depth = abs( $depth{$PE_cord{$r1.':'.$data[4]}} - $med_x ) / $mad_x;
      }
   elsif ($data[0] =~ /Y/)
      {
       $PE_depth = abs( $depth{$PE_cord{$r1.':'.$data[4]}} - $med_y ) / $mad_y;
      }
   else
      {
       $PE_depth = abs( $depth{$PE_cord{$r1.':'.$data[4]}} - $med ) / $mad;
      }

    next if ( ($PE_depth > 3) && ($FILTER>0) );
    
    if ($len > $CUT1*$size/100)
       {
        print OU "$line\n";
       }
    elsif ( ($len > $CUT2*$size/100) && ($te =~ "Alu") ) 
       {
        $seq1 = substr($data[12], 0, $data[8]);
        $seq2 = substr($data[12], $data[9], ($size-$data[9]));
        $pass = 0;
        $pass = 1 if ($seq1 && (($seq1 =~ /AAAAAA/) || ($seq1 =~ /TTTTTT/)) );
        $pass = 1 if ($seq2 && (($seq2 =~ /AAAAAA/) || ($seq2 =~ /TTTTTT/)) );

        next if (!$pass);
        print OU "$line\n";
       }
     else
       {
        print O2 "$line\n";
       }
    }
   
sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub mad
   {
    my @vals = @_;
    my $med = &median(@vals);
    my ($i, $mad, @diff);
    my $constant = 1.4826;
    for ($i=0; $i<=$#vals; $i++)
       {
        $diff[$i] = abs($vals[$i] - $med);
       }
    $mad = &median(@diff);
    return($mad*$constant);
   }


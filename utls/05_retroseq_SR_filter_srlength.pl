#! /usr/bin/perl -w

use strict;

### add depth filter for split reads support ###

my $in_file = $ARGV[1].'/'.$ARGV[2].'/'.$ARGV[0].'.sr.discover';
my $ou_file = $ARGV[1].'/'.$ARGV[2].'/'.$ARGV[0].'.sr.tabe.discover';
my $seg_file= $ARGV[1].'/'.$ARGV[2].'/seg.sr.scores';
my $FILTER  = $ARGV[3];
my $dep1_file = $ARGV[1].'/'.$ARGV[2].'/'.$ARGV[0].'.bed';
my $dep2_file = $ARGV[1].'/'.$ARGV[2].'/'.$ARGV[0].'.depth';
my ($r1, $SEG_CUT, $sim_ALU, $sim_L1, $sim_other, $CUT);

if ($FILTER == 0)
   {
    $SEG_CUT = 0;
    $sim_ALU = 85;
    $sim_L1  = 85;
    $sim_other = 85;
    $CUT = 0;
   }
elsif ($FILTER == 1)
   {
    $SEG_CUT = 1;
    $sim_ALU = 85;
    $sim_L1  = 85;
    $sim_other = 85;
    $CUT = 0;
   }
elsif ($FILTER == 2)
   {
    $SEG_CUT = 1;
    $sim_ALU = 92.67;
    $sim_L1  = 98;
    $sim_other = 95;
    $CUT = 200;
   }
(open (IN, "<$in_file")) || die "cannot open the in file\n";
(open (OU, ">$ou_file")) || die "cannot open the ou file\n";
(open (SEG,"<$seg_file")) || die "cannot open the seg file\n";
(open (DEP1, "<$dep1_file")) || die "cannot open the dep1 file\n";
(open (DEP2, "<$dep2_file")) || die "cannot open the dep2 file\n";

my ($line, @data);
my ($tag, $chr, $c1, $c2, $te, $ref, $count, %anno, $read, $cigar, $start, $end);
my (%seg, %uniq, %out, $i, $j, $k, $id);
my ($med, $mad, %depth, @dep_all, @dep_X, $med_x, $mad_x, @dep_Y, $med_y, $mad_y, %cord, @temp1, $depth);

while ($line = <DEP1>)
   {
    chomp($line);
    @data = split("\t", $line);
    @temp1 = split('__', $data[3]);
    $cord{$temp1[1]} = $data[0]."\t".$data[1]."\t".$data[2];
   }

$i = 0;
$j = 0;
$k = 0;
while ($line = <DEP2>)
   {
    chomp($line);
    @data = split("\t", $line);
    $depth{$data[0]."\t".$data[1]."\t".$data[2]} = $data[3] / ($data[2]-$data[1]);
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
$med   = &median(@dep_all);
$mad   = &mad(@dep_all);
$med_x = &median(@dep_X);
$mad_x = &mad(@dep_X);
$med_y = &median(@dep_Y);
$mad_y = &mad(@dep_Y);
print "$med\t$mad\t$med_x\t$mad_x\t$med_y\t$mad_y\n";

while ($line = <SEG>)
   {
    chomp($line);
    @data = split("\t", $line);
    next if ($data[1] < $SEG_CUT);
    $data[0] =~ /^\>(.+)\>(\w+)\>(\w+)\>(\w+)\>(\w+)\>(\w+)\>(\w+)\(/;
    $tag = $1;
    $i   = $2;
    $cigar = $3;
    $seg{$tag."\t".$i."\t".$cigar."\t".$4."\t".$5."\t".$6."\t".$7} = $data[1];
   }

while ($line = <IN>)
  {
   chomp($line);
   @data = split("\t", $line);
   ### filtering SR duplicate reads ###
   ### TE cord2 cigar 0 49 2196 2147 ###
   next if (!exists($seg{$data[4]."\t".$data[3]."\t".$data[6]."\t".$data[11]."\t".$data[12]."\t".$data[13]."\t".$data[14]}));
   $r1 = ($data[5] & 0x40) ? "r1" : "r2";
   $data[4] = $r1.':'.$data[4];
   $id = $data[8];

   if ($data[3] =~ /Alu/)
      {
       next if (($id <= $sim_ALU) && ($FILTER>0));
      }
   elsif ($data[3] =~ /L1/)
      {
       next if (($id <= $sim_L1) && ($FILTER>0));
      }
   else
      {
       next if (($id <= $sim_other) && ($FILTER>0));
      }

   next if ((!exists($cord{$data[4]})) || (!exists($depth{$cord{$data[4]}}))); 

   if ($data[0] =~ /X/)
      {
       $depth = abs( $depth{$cord{$data[4]}} - $med_x ) / $mad_x;
      }
   elsif ($data[0] =~ /Y/)
      {
       $depth = abs( $depth{$cord{$data[4]}} - $med_y ) / $mad_y;
      }
   else
      {
       $depth = abs( $depth{$cord{$data[4]}} - $med ) / $mad;
      }

#   print "$depth{$cord{$data[4]}}\t$depth\n" if ($data[4] =~ /75148895/);
   next if ( ($depth > 3) && ($FILTER>0) );
### checking orientation ... done ###
### checking transduction ... ok ###
### remove reads that rule out intact 3' end ###
   if ($data[6] =~ /(\d+)M\d+S$/)
      {
       $data[1] += $1;
       $data[2] += $1;
       $start = $data[13];
       $end   = $data[14];
       $te    = $data[3];
       if ( ( ($te =~ 'Alu' && $start < 270) || ($te =~ 'L1' && $start < 5900) ) && ($start > $end) )
          {
           $tag = 'BAD1';   ### no 3' ###
          }
       elsif ( ( ($te =~ 'Alu' && $end < 270) || ($te =~ 'L1' && $end < 5900) ) && ($start < $end) && ($data[12] < ($data[10] - 10) ) )
          {
           $tag = 'BAD2'; ### something inserted before the full length 3'end ###
          }
       elsif ( ($te =~ /L1/) && ( abs(($start + $end)/2-6064) < $CUT))
          {
           $tag = 'BAD3'; ### too short ###
          }
       else
          {
           $tag = 'OK';
          }
      }
    elsif ($data[6] =~ /^\d+S\d+M/)
      {
       $start = $data[14];
       $end   = $data[13];
       $te    = $data[3];
       if ( ( ($te =~ 'Alu' && $start < 270) || ($te =~ 'L1' && $start < 5900) ) && ($start > $end) )
          {
           $tag = 'BAD1';
          }
       elsif ( ( ($te =~ 'Alu' && $end < 270) || ($te =~ 'L1' && $end < 5900) ) && ($start < $end) && ($data[11] > 10 ) )
          {
           $tag = 'BAD2';
          }
       elsif ( ($te =~ /L1/) && ( abs(($start + $end)/2-6064) < $CUT))
          {
           $tag = 'BAD3'; ### too short ###
          }
       else
          {
           $tag = 'OK';
          }
       }
     else
       {
        $start = '';
        $end   = $tag = '';
       }
    $line = join("\t", @data);
    $tag = 'OK' if ($FILTER != 2);
    $out{$line."\t".$start."\t".$end."\t".$tag} = 1;
   }

if ($FILTER == 0)
   {
    for $line(keys %out)
       {
        print OU "$line\n";
       }
   }
else
   {
    %uniq = ();
    for $line(keys %out)
       {
        @data = split("\t", $line);
        $data[9] = abs($data[9]);
        $id = $data[0]."\t".$data[1]."\t".$data[3]."\t".$data[6]."\t".$data[9]."\t".$data[24]."\t".$data[25];
        if (!exists($uniq{$id}))
           {
            $uniq{$id} = $line; # if (!exists($uniq{$id}));
           }
        else
           {
            #print "$uniq{$id}\n$line\n\n"; # if ($uniq{$id} ne $line);
           }
       }
    for $id(keys %uniq)
       {
        print OU "$uniq{$id}\n";
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


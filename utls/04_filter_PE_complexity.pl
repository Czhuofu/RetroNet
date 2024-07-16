#! /usr/bin/perl -w

use strict;

####################################################################################
### filter complexity score                                                      ###
### unique read                                                                  ###
###    if the read is chimeric, save the read in another file                    ###
###    otherwise pick the read.te.chr.cord1 with the highest map%*map_length     ###
### checking the cords of the uniquely mapped end and the sequence of the TE end ###
####################################################################################

my $disc_file = $ARGV[1].'/'.$ARGV[2].'/'.$ARGV[0].'.alignfilter.discover';
my $seg_file  = $ARGV[1].'/'.$ARGV[2].'/seg.pe.scores';
my $out_file  = $ARGV[1].'/'.$ARGV[2].'/'.$ARGV[0].'.filter.discover';
my $ou2_file  = $ARGV[1].'/'.$ARGV[2].'/'.$ARGV[0].'.chimera.pe.reads';
my $FILTER    = $ARGV[3];

my ($te, $first, $seq, $read, $line, @data, %seg, %support, $dup, %duplicate);
my $CUT = 1;
my ($support, $qual, $dup_factor, %id, %all);
my (%chimera, @array, %top, %cand, %reads, @site1, @site2);
my ($c1, $s1, $s2, $site, $i, $j);

(open (DISC, "<$disc_file")) || die "cannot open the disc file\n";
(open (SEG , "<$seg_file ")) || die "cannot open the seg file\n";
(open (OUT , ">$out_file ")) || die "cannot open the out file\n";
(open (OU2 , ">$ou2_file ")) || die "cannot open the ou2 file\n";

while ($line = <SEG>)
  {
   chomp($line);
   @data = split("\t", $line);
   $data[0] =~ /^\>(.+)\(/;
   $seg{$1} = $data[1];
  }

while ($line = <DISC>)
  {
   chomp($line);
   @data = split("\t", $line);
   $read = $data[4].'###'.$data[3].'###'.$data[7].'###'.$data[8].'###'.$data[9];
   # print "$line\nchecking $read\n";
   if ($data[3] =~ /L1/)
      {
       $te = 'L1';
      }
   elsif ($data[3] =~ /Alu/)
      {
       $te = 'Alu';
      }
   elsif ($data[3] =~ /SVA/)
      {
       $te = 'SVA';
      }
   else
      {
       $te = $data[3];
      }

   ### one unique support read ###
   $support = $te."\t".$data[4]."\t".$data[0]."\t".$data[1];
   
   next if ( ($seg{$read} <= $CUT) && ($FILTER > 0));
   $data[3] = $te;
   $line = join("\t", @data);
   $qual = $data[7] * ($data[9] - $data[8]);

   $c1 = ($data[10] > $data[11]) ? $data[11] : $data[10]; # &min($data[10], $data[11]);
   $s1 = ($data[8] > $data[9]) ? $data[9] : $data[8]; #&min($data[8], $data[9]);
   $s2 = ($data[8] > $data[9]) ? $data[8] : $data[9];

   if (exists($id{$support}) && ($id{$support} < $qual) )
      { ### line is a better hit ###
       $all{$support} = $line;
       $id{$support}  = $qual;
      }
   elsif (!exists($id{$support}))
      {
       $all{$support} = $line;
       $id{$support}  = $qual;
      }

   $reads{$support}{$s1."\t".$s2} = 1;
   if (!exists($top{$support}))
      {
       $top{$support} = $c1;
      }
   else
      {
       if ( abs($c1 - $top{$support}) > 100)
          {
           $cand{$support} = 1;
          }
      }
   }

for $read(keys %cand)
   {
    $i = 0;
    @array = ();
    for $site(keys %{$reads{$read}})
       {
        $array[$i++] = $site;
       }
    for ($i=0; $i<$#array; $i++)
       {
        @site1 = split("\t", $array[$i]);

        for ($j=$i+1; $j<=$#array; $j++)
           {
            @site2 = split("\t", $array[$j]);
            if ( (abs($site1[1]-$site2[0]) < 10) || (abs($site2[1]-$site1[0]) < 10) )
               { ### it's a chimera ###
                print OU2 "$read\t$array[$i]\t$array[$j]\n";
                $chimera{$read} = 1;
               }
           }
        }
    }

for $support(keys %all)
   {
    print OUT "$all{$support}\n" if (!exists($chimera{$support}));
   }

sub reverse_complement_IUPAC {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}


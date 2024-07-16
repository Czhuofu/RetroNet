#! /usr/bin/perl -w
#
use strict;
### need to consider TE names in locating PCR duplicates ###
### fixed a bug for changes made in step 06 ###
my $sub = $ARGV[0];
my $retro = $ARGV[1];
my $te = $ARGV[2];
my $path = $ARGV[3];
my $FILTER = $ARGV[4];
my $call_file = $path.'/'.$retro.'/'.$te.'/'.$sub.'.'.$te.'.PE.calls';
my $site_file = $path.'/'.$retro.'/'.$te.'/'.$sub.'.'.$te.'.novel.sites';
my $disc_file = $path.'/'.$retro.'/'.$sub.'.alignfilter.discover';
my $dup_file  = $path.'/'.$retro.'/'.$te.'/'.$sub.'.'.$te.'.dup.reads';
my $clust_file= $path.'/'.$retro.'/'.$te.'/'.$sub.'.'.$te.'.dup.cluster.txt'; 
my $new_file  = $path.'/'.$retro.'/'.$te.'/'.$sub.'.'.$te.'.PE.nodup.calls';
my $new_site  = $path.'/'.$retro.'/'.$te.'/'.$sub.'.'.$te.'.nodup.sites';
my ($line, @data, @temp, %dup, %num, %site, %count, %pair);
my (@pe, @call, %label, %remove, @sites);
my ($id, $i, $j, $k, $same1, $same2, $same3, $first, $count, $subte, $tetag);
my %TEtag = (
   'LINE' => 'L1',
   'ALU'  => 'Alu',
   'SVA'  => 'SVA',
   'HERVK' => 'HERV', 
   'HERVH' => 'HERV', 
  );
$tetag = $TEtag{$te};

if ($FILTER==0)
   {
    system ("cp $call_file $new_file");
    exit;
   }

(open (CALL, "<$call_file")) || die "cannot open the call file\n";
(open (SITE, "<$site_file")) || die "cannot open the site file\n";
(open (DISC, "<$disc_file")) || die "cannot open the disc file\n";
(open (DUP , ">$dup_file ")) || die "cannot open the dup file\n";
(open (NEW , ">$new_file ")) || die "cannot open the dup file\n";
(open (NSITE, ">$new_site")) || die "cannot open the new site file\n";
(open (CLUST,">$clust_file")) || die "cannot open the clust file\n";

%num = ();
@call = ();
$j   = 0;
while ($line = <CALL>)
   {
    chomp($line);
    $call[$j++] = $line;
    @data = split("\t", $line);
    next if ($data[3] < 2);
    $id   = $data[0]."\t".$data[1]."\t".$data[2];
    for ($i=4; $i<=$#data; $i++)
       {
        $num{$data[$i]} = $id;
       }
   }
$i = scalar keys %num;
print "read in all reads in num $i\n";

$i = 0;
while ($line = <SITE>)
   {
    chomp($line);
    $i ++;
    $sites[$i] = $line;
    if (exists($num{$i}))
       { ### this read is supporting an insertion ###
        @data = split("\t", $line);
        $id = $data[4]."\t".$data[14];  ### read ID and flag ###
        $label{$id} = $i; ### read name for support read $i  ### 
       }
   }

my %all = ();
while ($line = <DISC>)
   {
    chomp($line);
    @data = split("\t", $line);
    $subte = $data[3];
    next if ($subte !~ /$tetag/);
    $id = $data[4]."\t".$data[14];
    if ( exists($label{$id}) && (!exists($all{$data[4]."\t".$data[14]."\t".$subte})) )
       {
        $i = $label{$id};
        $count{$num{$i}}{$subte} = 0 if (!exists($count{$num{$i}}{$subte})); ### number of PE support reads for chr.c1.c2 insertion ###
        $count = $count{$num{$i}}{$subte};
        if ($data[4] =~ /\_/)
           {
            @temp = split('_', $data[4]);
            $site{$num{$i}}{$subte}[$count][0] = $temp[0];  ### library ID for the support read $i ###
           }
        else
           {
            $site{$num{$i}}{$subte}[$count][0] = 'lib';
           }
        $site{$num{$i}}{$subte}[$count][1] = $data[1];     ### hgcord1 ###  
        $site{$num{$i}}{$subte}[$count][2] = $data[2];     ### hgcord2 ###
        if ($data[15] eq '+')
           {
            $site{$num{$i}}{$subte}[$count][3] = $data[10];  ### tecord1 ###
            $site{$num{$i}}{$subte}[$count][4] = $data[11];  ### tecord2 ###
           }
        else
           {
            $site{$num{$i}}{$subte}[$count][4] = $data[10];  ### tecord2 ###
            $site{$num{$i}}{$subte}[$count][3] = $data[11];  ### tecord1 ###
           }
        $site{$num{$i}}{$subte}[$count][5] = $data[4]."\t".$data[14]; ### read name ###
        $site{$num{$i}}{$subte}[$count][6] = $data[16]; ### CIGAR ###
        $all{$data[4]."\t".$data[14]."\t".$subte} = 1;
        $count{$num{$i}}{$subte} ++;
       }
   }

for $id(keys %site)
   {
    for $subte(keys %{$site{$id}})
       {
        for ($i=0; $i<($count{$id}{$subte}-1); $i++)
           {
            for ($j=$i+1; $j<$count{$id}{$subte}; $j++)
               {
                $same1 = $same2 = 0;
                next if ($site{$id}{$subte}[$i][0] ne $site{$id}{$subte}[$j][0]);
                $same3 = &parse_anchor($site{$id}{$subte}[$i][6], $site{$id}{$subte}[$j][6]);
                $same1 ++ if (($site{$id}{$subte}[$i][2] eq $site{$id}{$subte}[$j][2]) || ($site{$id}{$subte}[$i][1] eq $site{$id}{$subte}[$j][1]));
                $same2 ++ if (($site{$id}{$subte}[$i][4] eq $site{$id}{$subte}[$j][4]) || ($site{$id}{$subte}[$i][3] eq $site{$id}{$subte}[$j][3]));

#if ( ($site{$id}{$subte}[$i][5] =~ /$ARGV[5]/) && ($site{$id}{$subte}[$j][5] =~ /$ARGV[6]/) )
#   {
#    print "$same1\t$same2\t$same3\t$site{$id}{$subte}[$i][6]\t$site{$id}{$subte}[$j][6]\n";
#   }
#elsif ( ($site{$id}{$subte}[$i][5] =~ /$ARGV[6]/) && ($site{$id}{$subte}[$j][5] =~ /$ARGV[5]/) )
#   {
#    print "$same1\t$same2\t$same3\t$site{$id}{$subte}[$i][6]\t$site{$id}{$subte}[$j][6]\n";
#   }
               if ($same1 && $same2 && $same3)
                  {
                   $dup{$id}{$site{$id}{$subte}[$i][5]} = 1;   ### read name is an duplicate ###
                   $dup{$id}{$site{$id}{$subte}[$j][5]} = 1;
                   $pair{$id}{$site{$id}{$subte}[$i][5]}{$site{$id}{$subte}[$j][5]} = 1; ### pairs of all duplcates ###
                   print DUP "$id\t$subte\t$site{$id}{$subte}[$i][0]\t$site{$id}{$subte}[$i][5]\t$site{$id}{$subte}[$i][1]\t$site{$id}{$subte}[$i][2]\t$site{$id}{$subte}[$i][3]\t$site{$id}{$subte}[$i][4]\t";
                   print DUP "$site{$id}{$subte}[$j][5]\t$site{$id}{$subte}[$j][1]\t$site{$id}{$subte}[$j][2]\t$site{$id}{$subte}[$j][3]\t$site{$id}{$subte}[$j][4]\n";
                  }
               }
            }
        }
    }
my @all = ();
my %assign = ();
my ($name);
#my $cluster_i = 0;

for $id(keys %dup)
   {
    @all = ();
    %assign = ();
    $i = 0;
    for $name(keys %{$dup{$id}})
       { ### all reads with at least one duplicate ###
        $all[$i++] = $name;
       }

    $i = 0;
    for ($j=0; $j<$#all; $j++)
      {
       if (!exists($assign{$all[$j]}))
         {
          $assign{$all[$j]} = $i; ### assign dup reads to clusters ###
          $i ++;
         }
       for ($k=$j+1; $k<=$#all; $k++)
         {
          if (exists($pair{$id}{$all[$k]}{$all[$j]}) || exists($pair{$id}{$all[$j]}{$all[$k]}) )
            {
             $assign{$all[$k]} = $assign{$all[$j]} if (!exists($assign{$all[$k]}));
            }
         }
      }

   for ($j=0; $j<$i; $j++)
      { ### check each cluster ###
       print CLUST "$id\t$j";
       $first = 1;
       for ($k=0; $k<=$#all; $k++)
          {
           if ($assign{$all[$k]} eq $j)
              {
               print CLUST "\t$label{$all[$k]}\|$all[$k]";
               if ($first == 0)
                  {
                   $remove{$label{$all[$k]}} = 1; ### remove all reads except for the first ###
                  }
               $first = 0;
              }
          }
       print CLUST "\n";
      }
   }

### print new sites ###
for ($j=1; $j<=$#sites; $j++)
   {
    print NSITE "$sites[$j]\n" if (!exists($remove{$j}));
   }

### print new output ###
for ($j=0; $j<=$#call; $j++)
   {
    @data = split("\t", $call[$j]);
    if ($data[3] < 2)
       {
        print NEW "$data[0]\t$data[1]\t$data[2]\t$te\t$data[3]\t$data[4]\n";
       }
    else
       {
        print NEW "$data[0]\t$data[1]\t$data[2]\t$te";
        $count=0;
        for ($i=4; $i<=$#data; $i++)
           {
            $count ++ if (!exists($remove{$data[$i]}));
           }
        print NEW "\t$count";
        for ($i=4; $i<=$#data; $i++)
           {
            print NEW "\t$data[$i]" if (!exists($remove{$data[$i]}));
           }
        print NEW "\n";
        }
    }

sub parse_anchor
    {
     my $cigar1 = shift;
     my $cigar2 = shift;
     my ($c1, $c2, $c3, $c4, $flag);
     $flag = 1;
     if (($cigar1 !~ /S/) && ($cigar2 !~ /S/))
       {
        $flag = 1;
       }
     elsif (($cigar1 =~ /S/) && ($cigar2 !~ /S/))
       {
        $flag = 0;
       }
     elsif (($cigar1 !~ /S/) && ($cigar2 =~ /S/))
       {
        $flag = 0;
       }
     elsif ($cigar1 =~ /(\d+)M(\d+)S/)
       {
        $c1 = $1;
        $c2 = $2;
        if ($cigar2 =~ /(\d+)M(\d+)S/)
          {
           $c3 = $1;
           $c4 = $2;
           $flag = 0 if ( ($c1 != $c3) && ($c2 != $c4) );
          }
        else
          {
           $flag = 0;
          }
       }
     elsif ($cigar1 =~ /(\d+)S(\d+)M/)
       {
        $c1 = $1;
        $c2 = $2;
        if ($cigar2 =~ /(\d+)S(\d+)M/)
          {
           $c3 = $1;
           $c4 = $2;
           $flag = 0 if ( ($c1 != $c3) && ($c2 != $c4) );
          }
        else
          {
           $flag = 0;
          }
       }
     return($flag);
    }

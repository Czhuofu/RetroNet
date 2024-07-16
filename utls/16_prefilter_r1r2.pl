#!/usr/bin/perl -w

use strict;

### keep one read for both r1 and r2 for each MEI ###
#### SR read over PE read                         ###
my $sub  = $ARGV[0];
my $path = $ARGV[1];
my $retro= $ARGV[2];
my $TE   = $ARGV[3];
my $in_file = $path.'/'.$retro.'/'.$TE.'/'.$sub.'.'.$TE.'.noref.calls';
my $ou_file = $path.'/'.$retro.'/'.$TE.'/'.$sub.'.'.$TE.'.novel.calls';
my $all_PE  = $path.'/'.$retro.'/'.$TE.'/'.$sub.'.'.$TE.'.novel.sites';;

my %tag = (
   'LINE' => 'L1' ,
   'ALU'  => 'Alu',
  );
my $tetag = $tag{$TE};
my ($line, @data);
my ($i, $j, $r1, $pe_num, $sr_num, $total_num, $key, $sr, $sr_read);
my (@allpe, %pe_reads, %sr_reads, @temp, @temp2, @nums, @dat2, @sr);

(open (IN, "<$in_file")) || die "cannot open $in_file\n";
(open (OU, ">$ou_file")) || die "cannot open $ou_file\n";
(open (APE, "<$all_PE")) || die "cannot open $all_PE\n";

$i = 0;
while ($line = <APE>)
   {
    chomp($line);
    $i ++;
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $allpe[$i][0] = $data[4];
    if ($data[4] =~ /^(lib\d+)\_(.+)/)
       {
        $allpe[$i][1] = $1.'_'. $r1.':'.$2;
       }
    else
       {
        $allpe[$i][1] = $r1.':'.$data[4];
       }
   }

while ($line = <IN>)
   {
    chomp($line);
    @data = split("\t", $line);
    $pe_num = $sr_num = $total_num = 0;
    %pe_reads = %sr_reads = ();

    if ( ($data[7] =~ /pe\,/) && ($data[7] =~ /sr\,/) )
       {
        @temp = split(';', $data[7]);
        for ($j=0; $j<=$#temp; $j++)
            {
             if ($temp[$j] =~ /pe/)
                {
                 @nums = split(',', $temp[$j]);
                 for ($i=4; $i<=$#nums; $i++)
                    {
                     if ($nums[$i])
                       {
                        @dat2 = split('_', $nums[$i]);
                        $pe_reads{$dat2[0]} = 1;
                       }
                    }
                }
             elsif ($temp[$j] =~ /sr/)
                {
                 @nums = split(',', $temp[$j]);
                 @sr   = split(/\|/, $nums[2]);
                 foreach $sr(@sr)
                    {
                     ($r1, @temp2) = split(':', $sr);
                     $sr_read = join(':', @temp2);
                     $sr_reads{$sr_read} = $sr;
                    }
                }
            }
       }
    elsif ($data[7] =~ /pe\,/)
       {
        @nums = split(',', $data[7]);
        for ($i=4; $i<=$#nums; $i++)
           {
            if ($nums[$i] && ($nums[$i] =~ /\_/))
              {
               @dat2 = split('_', $nums[$i]);
               $pe_reads{$dat2[0]} = 1;
              }
          }
       }
    else
       {
        @sr = split(/[\,\|\;]/, $data[7]);
        foreach $sr(@sr)
          {
           next if (($sr eq 'sr') || ($sr =~ /^\d+$/) || ($sr !~ /\w/));
           ($r1, @temp2) = split(':', $sr);
           $sr_read = join(':', @temp2);
           $sr_reads{$sr_read} = $sr;
          }
       }

   $sr_num = scalar keys %sr_reads;
   for $key(keys %pe_reads)
      {
       if ((scalar keys %sr_reads > 0 ) && exists($sr_reads{$allpe[$key][0]}))
          {
           delete $pe_reads{$key};
          }
      }
   $pe_num = scalar keys %pe_reads;
   $total_num = $pe_num + $sr_num;
   print OU "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$total_num\t$sr_num\t$pe_num\t";
   if ($pe_num > 0)
      {
       print OU "pe";
       for $key(keys %pe_reads)
          {
           print OU "\,$allpe[$key][1]";
          }
       if ($sr_num > 0)
          {
           print OU "\;sr";
           for $sr_read(keys %sr_reads)
              {
               print OU "\,$sr_reads{$sr_read}";
              }
           }
       }
    elsif ($sr_num > 0)
       {
        print OU "sr";
        for $sr_read(keys %sr_reads)
           {
            print OU "\,$sr_reads{$sr_read}";
           }
        }
    print OU "\n";
  }



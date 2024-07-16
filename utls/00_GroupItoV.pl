#! /usr/bin/perl -w

use strict;

### split PE reads to PE and SR ###
### fix $split, no longer add $M ###
### fix the flag for SR2 ###
my $masterpath = $ARGV[0];
my $subject = $ARGV[1];
my $workfolder = $ARGV[2];
my $pe_file = $workfolder.'/'.$subject.'.discover';
my $sr_file = $workfolder.'/'.$subject.'.sr.discover';
my $sr2_file = $workfolder.'/'.$subject.'.sr2.discover';
my $old_pe  = $workfolder.'/'.$subject.'.old.discover';
my $old_sr  = $workfolder.'/'.$subject.'.sr.old.discover';
my $clipFasta = $workfolder.'/'.$subject.'.sr.remap.fa';

### backing up the PE and SR supporting reads ###
system("mv $pe_file $old_pe");
system("cp $sr_file $old_sr");

my ($gap, $c1, $c2, $c3, $c4, $line, @data);
my ($te, $read, $flag, $newflag, $SA, $id, $split, $seq, $map, $distance, $i, $strand);
my ($lastLine, $input, %anno, @s, %seq, @old, %move, @temp);
my $count = 0;
my $OVER  = 10;

(open (OL, "<$old_pe")) || die "cannot open the old file\n";
(open (PE, ">$pe_file")) || die "cannot open the pe file\n";
(open (SR, ">>$sr_file")) || die "cannot open the sr file\n";
(open (SR2, ">$sr2_file")) || die "cannot open the sr2 file\n";
(open (FA, ">$clipFasta")) || die "cannot open the clipfasta\n";

while($line = <OL>)
   {
    chomp($line);
    @data = split("\t", $line);
    $old[$count++] = $line;
    ($flag, $SA) = &regroup($data[0], $data[1], $data[2], $data[8], $data[9], $data[13], $data[14], $data[19], $data[21]);
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

    if ($flag)
       {
        # print "$flag $SA\n";
        $id = $te.'###'.$data[4].'###'.$data[14];
        ###  chr1,2683515,+,75M76S,6,3 ###
        @temp = split(',', $SA);
        $move{$data[4]} = 1;
        if ($temp[2] eq '+')
           { ### positive strand ###
            if ($data[14] < 128)
               { ### insert in on read 1
                $newflag = 1 + 2 + 32 + 64;  ### 99 read1 + ###
               }
            else
               {### insert in on read 2
                $newflag = 1 + 2 + 32 + 128; ### 163 read2 + ###
               }
            if ($data[15] eq '+')
               {
                $read = $data[12];
                $c1 = $data[8];
                $c2 = $data[9];
                $c3 = $data[10];
                $c4 = $data[11];
                $strand = $data[6];
               }
            else
               { ### SA is on another strand ###
                $c1 = length($data[12]) - $data[9];
                $c2 = length($data[12]) - $data[8];
                $c3 = $data[11];
                $c4 = $data[10];
                $read = &rc($data[12]);
                $strand = ($data[6] eq '+') ? '-' : '+';
               } 
           }
        else
           { ### negative strand ###
            if ($data[14] < 128)
               { ### insert in on read 1
                $newflag = 1 + 2 + 16 + 64;   ### 83 read1 - ###
               }
            else
               {### insert in on read 2
                $newflag = 1 + 2 + 16 + 128;  ### 147 read2 - ###
               }
            if ($data[15] eq '-')
               {
                $read = $data[12];
                $c1 = $data[8];
                $c2 = $data[9];
                $c3 = $data[10];
                $c4 = $data[11];
                $strand = $data[6]; # eq '+') ? '-' : '+';
               }
            else
               {
                $c1 = length($data[12]) - $data[9];
                $c2 = length($data[12]) - $data[8];
                $c3 = $data[11];
                $c4 = $data[10];
                $read = &rc($data[12]);
                $strand = ($data[6] eq '+') ? '-' : '+';
               }
           }
    
        if ($temp[3] =~ /^(\d+)M(\d+)S$/)
           {
            $split = $temp[1]; # + $1;
            $seq = substr($read, $1, $2);
           }
        elsif ($temp[3] =~ /^(\d+)S(\d+)M$/)
           {
            $split = $temp[1];
            $seq = substr($read, 0, $1);
           }
        
        $distance = $temp[1] - $data[1];
        $map = $data[9] - $data[8];

        ### $seq uniquely correlates with TE_Read_Flag ### 
        $anno{$id} = "$temp[0]\t$split\t$split\t$data[3]\t$data[4]\t$newflag\t$temp[3]\t$strand\t$data[7]\t$distance\t$map\t$c1\t$c2\t$c3\t$c4\t$seq\t$seq\t$read\t$data[16]\t$data[17]\t$data[18]\t$data[19]\t$data[22]\t$data[23]\t$data[24]\t$data[1]"; 
        #print "exonerate ID $id\n$anno{$id}\n" if ($id =~ /1:2101:11459:151891/);
        $seq{$id."\t\t".$seq} = 1;
       }
   }

for $line(sort keys %seq)
   {
    ($id, $seq) = split("\t\t", $line);
    print FA ">$id\n$seq\n";
   }
close(FA);
### exonerate ###
my $exo_file = $workfolder.'/exonerate.SR.out';
my $refsFasta = $masterpath.'/refTE/allrefs.fasta';
system("env TMPDIR=$workfolder exonerate --model affine:local --bestn 1 --ryo \"INFO: %qi %qal %pi %tS %ti %ql %qab %qae %tab %tae %tas %qas %qs\\n\" $clipFasta $refsFasta > $exo_file"); 
my %uniq = ();
(open (EXO, "<$exo_file")) || die "cannot open the exo_file\n";
while ($line = <EXO>)
   {
    chomp($line);
    $lastLine = $line;
    if ($line =~ /vulgar/)
       {
        $gap = 0;
        while ($line =~ /G (\d+) (\d+)/g)
           {
            $gap += $1;
            $gap += $2;
           }
        }
    if ($line =~ /^INFO/)
       {
        $input = $line;
        while ( ($line = <EXO>) && ($line =~ /^[ACTGN actgn]+/) && ($line !~ /^C4/) )
           {
            chomp($line);
            $input .= $line;
           }
        $input = $input." ".$gap;
        @data = split(" ", $input);
        $id = $data[1];
        if (exists($anno{$id}))
           {
            #print "exonerate ID $id\n" if ($id =~ /1:2101:11459:151891/);
            @s = split("\t",$anno{$id});
            $line = "$s[0]\t$s[1]\t$s[2]\t$data[5]\t$s[4]\t$s[5]\t$s[6]\t$data[4]\t$data[3]\t$s[9]\t$data[6]\t$data[7]\t$data[8]\t$data[9]\t$data[10]\t$data[11]\t$data[12]\t$data[13]\t$s[18]\t$s[19]\t$s[20]\t$s[21]\t$s[22]\t$s[23]\t$s[24]\t$s[25]\t$data[14]";
            if (!exists($uniq{$line}))
               {
                print SR2 "$line\n";
                print SR "$s[0]\t$s[1]\t$s[2]\t$data[5]\t$s[4]\t$s[5]\t$s[6]\t$s[7]\t$data[3]\t$s[9]\t$data[6]\t$data[7]\t$data[8]\t$data[9]\t$data[10]\t$data[11]\t$data[12]\t$data[13]\t$s[21]\t$s[18]\t$s[20]\t$s[22]\t$s[23]\t$s[24]\t$s[25]\t$data[14]\n";
                $uniq{$line} = 1;
               }
           }
        }
    }

for ($i=0; $i<=$#old; $i++)
   { ### reads cannot be both PE and SR supports ###
    @data = split("\t", $old[$i]);
    $id = $data[4];
    if (!exists($move{$id}))
       {### reads in %move are now considered properly paired, and will be removed from PE candidates ###
        print PE "$old[$i]\n";
       }
   }

sub regroup
   {
    my $chr = shift;
    my $c1  = shift;
    my $c2  = shift;
    my $te1 = shift;
    my $te2 = shift;
    my $cigar= shift;
    my $flag = shift;
    my $mate = shift;
    my $SA = shift;
    my $INSERT = 3000;
    my $overhang = 20;

    my ($tag, @tags, $mate_strand, @temp);
    $mate_strand = ($mate & 0x10) ? '-' : '+';
    if ($SA eq 'SA_NA')
       {
        return(0, '');
       }
    else
       { 
        ### get all SA alignments ###
        ###  chr1,2683515,+,30S45M76S,6,3;chrUn_KN707905v1_decoy,758,-,75S76M,34,3; ###
        @tags = split(';', $SA);

        foreach $tag(@tags)
           {
            ### properly paired ? ###
            @temp = split(',', $tag);
            ### same chr, close coordinates and different strands ###
            #print "$chr $temp[0] $temp[1] $c1 $temp[2] $mate_strand\n";
            next if ( ($chr ne $temp[0]) || (abs($temp[1] - $c1)> $INSERT) || ($temp[2] eq $mate_strand) );
            ### Cigar is chimeric to SA and to TE cords ###
            if ( ($temp[3]=~tr/S/S/) == 1) #&& &chimer($cigar, $temp[3], $te1, $te2) )
               {
                if ($temp[3] =~ /^(\d+)M(\d+)S$/)
                   {
                    return(1, $tag) if ($2 > $overhang);
                   }
                elsif ($temp[3] =~ /^(\d+)S(\d+)M$/)
                   {
                    return(1, $tag) if ($1 > $overhang);
                   }
               }
           }
       }
    return(0, '');
   }

sub rc 
  { ### reverse complement ###
   my $seq = shift;
   my %RCbase = (
      'A' => 'T',
      'T' => 'A',
      'G' => 'C',
      'C' => 'G',
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
   return($seqout);
  }


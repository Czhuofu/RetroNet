#!/usr/bin/perl -w

use strict;

my $BAMFLAGS =
{
    'paired_tech'    => 0x0001,
    'read_paired'    => 0x0002,
    'unmapped'       => 0x0004,
    'mate_unmapped'  => 0x0008,
    'reverse_strand' => 0x0010,
    'mate_reverse'   => 0x0020,
    '1st_in_pair'    => 0x0040,
    '2nd_in_pair'    => 0x0080,
    'not_primary'    => 0x0100,
    'failed_qc'      => 0x0200,
    'duplicate'      => 0x0400,
};

my $support=&isSupportingClusterRead(163,152,25,20,20,"57S44M","","","");

print "support=$support\n";

sub isSupportingClusterRead
{
    die qq[Incorrect number of parameters: ].scalar(@_) unless @_ == 9;

    my $flag = shift;
    my $insert = shift;
    my $mapQ = shift;
    my $minQual = shift;
    my $minSoftClip = shift;
    my $cigar = shift;
    my $refpos = shift; #if the breakpiont position + read start is provided, then check if the reads soft clip OVER the breakpoint
    my $readpos = shift;
    my $readlength = shift;
print "testing flag=$flag\ninsert=$insert\nmapQ=$mapQ\nminQual=$minQual\nSoftClip=$minSoftClip\ncigar=$cigar\nrefpos=$refpos\nreadpos=$readpos\nlength=$readlength\n";
    #            read is not a duplicate        map quality is >= minimum
    if( ! ( $flag & $$BAMFLAGS{'duplicate'} ) && $mapQ >= $minQual )
    {
print "not a duplicate and has good mapQ\n";
        #           read is mapped                       mate is unmapped
        if( ! ( $flag & $$BAMFLAGS{'unmapped'} ) && ( $flag & $$BAMFLAGS{'mate_unmapped'} ) )
        {
print "read is mapped, mate is unmapped\n";
            return 1;
        }
        #                               read mapped                             mate mapped                            ins size sensible          has single soft clip in cigar bigger than minimum clip size
        elsif( defined( $minSoftClip ) && $minSoftClip > 0 && ! ( $flag & $$BAMFLAGS{'unmapped'} ) && ! ( $flag & $$BAMFLAGS{'mate_unmapped'} ) && ( $flag & $$BAMFLAGS{'read_paired'} ) && abs( $insert ) < 3000 && ($cigar=~tr/S/S/) == 1 && $cigar =~ /(\d+)(S)/ && $1 > $minSoftClip )
        {
print "read is mapped, has a big clip\n";
            #if breakpoint information is provided - check the soft clip is consistent with the breakpoint
            if( defined( $refpos ) && defined( $readpos ) )
            {
print "soft clip consistent with the breakpoint\n";
                #                               forward aligned and 3' of break         cigar ends with match
                if( ( !( $flag & $$BAMFLAGS{'reverse_strand'} ) && $readpos > $refpos && $cigar =~ /M$/ && abs($refpos-$readpos)<$readlength ) || (( $flag & $$BAMFLAGS{'reverse_strand'} ) && $readpos < $refpos && $cigar =~ /S$/ && abs($refpos-$readpos)<$readlength ) )
                {
print "forward aligned and 3' of break         cigar ends with match\n";
                    return 3;
                }
                else{return 0;}
            }
            return 3;
        }
        #            read is mapped                         mate is mapped                                  not paired correctly                has large enough deviation from expected ins size (short libs assumed)
        elsif( ! ( $flag & $$BAMFLAGS{'unmapped'} ) && ! ( $flag & $$BAMFLAGS{'mate_unmapped'} ) && ! ( $flag & $$BAMFLAGS{'read_paired'} ) && ( abs( $insert ) > 30000 || $insert == 0 ) )
        {
print "read is mapped, mate is mapped, but not correctly paired\n";
            return 2;
        }
    }
    return 0;
}



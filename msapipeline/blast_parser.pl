#!/usr/bin/perl
#
#
# BLAST Parser  -  Script for parsing the BLAST output
#
# Version 1.1.5 (June 27, 2012)
#
# Copyright (c) 2011-2012 Kirill Kryukov
#
# This software is provided 'as-is', without any express or implied
# warranty. In no event will the authors be held liable for any damages
# arising from the use of this software.
#
# Permission is granted to anyone to use this software for any purpose,
# including commercial applications, and to alter it and redistribute it
# freely, subject to the following restrictions:
#
# 1. The origin of this software must not be misrepresented; you must not
#    claim that you wrote the original software. If you use this software
#    in a product, an acknowledgment in the product documentation would be
#     appreciated but is not required.
# 2. Altered source versions must be plainly marked as such, and must not be
#    misrepresented as being the original software.
# 3. This notice may not be removed or altered from any source distribution.
#


#
# This script is a BLAST parser - it will process the BLAST output and convert
# it into a compact form where each hit is represented by a single line.
#
# Usage:
#   perl blast_parser.pl <blastout.txt >parsed.txt
#
# You can also use it to parse the BLAST output on the fly as it is generated,
# saving a lot of disk space (although losing all alignments).
#

#> VIR1119
#Length=3215

# Score =  675 bits (365),  Expect = 0.0
# Identities = 396/411 (96%), Gaps = 2/411 (0%)
# Strand=Plus/Minus

use strict;
my $PRINTALL=0;#$ARGV[0];

my ($qname,$sname) = ('','');
my ($qchanged,$schanged,$hchanged) = (0,0,0);
my ($qnaming,$snaming,$qlen,$slen) = (0,0,0,0);
my ($qstart,$qend,$sstart,$send) = (-1,-1,-1,-1);
my ($bits,$expect,$length,$ident,$positives,$gaps,$frame,$HStrand,$TStrand) = (-1,'',-1,-1,-1,-1,'',"","");
my ($n_queries,$n_hits) = (0,0);

print "#\@SEQUENCE\n#Sequence Length\n";
print "#>genome\tgenome_length\tident\tgaps\tbits\texpect\tqstart\tqend\tsstart\tsend\n";
while (<>)
{
    chomp;

    if ($qnaming)
    {
        if ($_ =~ /^\s*Length=(\d+)$/ or $_ =~ /^\s*\((\d+) letters\)$/)
        {
            $qlen = $1;
            $qnaming = 0;
            #print STDERR "Query = $qname\n";
            #print STDERR ">";
        }
        else { $qname .= " $_"; }
        next;
    }

    if ($snaming)
    {
        if (/^\s*Length\s*=\s*(\d+)$/)
        {
            $slen = $1;
            $snaming = 0;
            #print STDERR "    Subject = $sname\n";
            #print STDERR ":";
        }
        else { s/^\s+//; $sname .= " $_"; }
        next;
    }

    if (/^Query[\:\s]\s(\d+)\s.*\s(\d+)$/)
    {
        if ($qstart < 0 or $1 < $qstart) { $qstart = $1; }
        if ($qend < 0 or $1 > $qend) { $qend = $1; }
        if ($2 < $qstart) { $qstart = $2; }
        if ($2 > $qend) { $qend = $2; }
        next;
    }

    if (/^Sbjct[\:\s]\s(\d+)\s.*\s(\d+)$/)
    {
#if ($sstart < 0 or $1 < $sstart) { $sstart = $1; }
#        if ($send < 0 or $1 > $send) { $send = $1; }
#        if ($2 < $sstart) { $sstart = $2; }
#        if ($2 > $send) { $send = $2; }

# Keep blast map direction intact..
	if($sstart<0) {$sstart=$1;}
	$send=$2;

        next;
    }

    if (/^\sScore/)
    {
        if ($snaming) { print STDERR "Error: can't complete subject name: $sname\n"; next; }
        if ($qnaming) { print STDERR "Error: can't complete query name: $qname\n"; next; }
        &flush_hit();
        ($qstart,$qend,$sstart,$send) = (-1,-1,-1,-1);
        ($bits,$expect,$length,$ident,$positives,$gaps,$frame) = (-1,'',-1,-1,-1,-1,'');

        if (/^\sScore\s+=\s+([\d\+e\.]+)\s+bits/) { $bits = int($1); } else { print STDERR "Error: Can't parse score: \"$_\"\n"; }
        if (/\s+Expect\(?\d*\)?\s+=\s+([e\-0-9\.]+)/) { $expect = $1; } else { print STDERR "Error: Can't parse expectation: \"$_\"\n"; }
        $hchanged = 1;
        next;
    }

    if (/^\sIdentities/)
    {
        ($length,$ident,$positives,$gaps,$frame) = (-1,-1,-1,-1,'');
        if (/^\sIdentities\s+=\s+(\d+)\/(\d+)\s+\(\d+\%\)/)
        {
            ($ident,$length) = ($1,$2);
        }
        if (/\s+Positives\s+=\s+(\d+)\/\d+\s+\(\d+\%\)/)
        {
            ($positives) = ($1);
        }
        if (/Gaps\s+=\s+(\d+)\//)
        {
            $gaps = $1;
        }
        else
        {
            $gaps = 0;
        }
        if ($length < 0)
        {
            print STDERR "Error: unexpected BLAST ' Identities = ' line format: $_\n";
        }
        next;
    }

    if (/\s[Ff]rame\s=\s(\S+)$/)
    {
        $frame = $1;
        if ($frame !~ /^(-?\+?\d|-?\+?\d\/-?\+?\d)$/) { print STDERR "Error: Unexpected frame: \"$frame\"\n"; }
        next;
    }

    if (/\sStrand=/)
    {
        my $Line = $_;
	{
		chomp $Line;
		my @PM=split(/[=\/]/,$Line);
		if($PM[1] eq "Plus")
		{
			$HStrand="+";
		}
		else
		{
			if($PM[1] ne "Minus")
			{
				die "Parse error..";
			}
			$HStrand="-";
		}
		if($PM[2] eq "Plus")
		{
			$TStrand="+";
		}
		else
		{
			if($PM[2] ne "Minus")
			{
				die "Parse error..";
			}
			$TStrand="-";
		}
	}
        next;
    }

    if (/^Query=\s+(.+)$/)
    {
        my $new_qname = $1;
        &flush_hit();
        $qname = $new_qname;
        $qnaming = 1;
        $qchanged = 1;
        $n_queries++;
        next;
    }

    if (/^>(.+)$/)
    {
        my $new_sname = $1;
        &flush_hit();
        $sname = $new_sname;
        $snaming = 1;
        $schanged = 1;
        next;
    }
}

&flush_hit();
print STDERR "OK, processed $n_queries queries, $n_hits hits\n";

sub flush_hit
{
    unless ($hchanged) { return; }
    if ($qstart<0 || $qend<0 || $sstart<0 || $send<0)
    {
        print STDERR "Error: no coordinates!\n";
        $hchanged = 0;
        return;
    }
    if ($bits<0 || $expect eq '' || $length<0 || $ident<0)
    {
        print STDERR "Error: parameters missing!\n";
        $hchanged = 0;
        return;
    }
    if ($frame eq '') { $frame = '-'; }

    #if ($frame eq '')
    #{
    #    print STDERR "Error: Frame/Strand is missing!\n";
    #    $hchanged = 0;
    #    return;
    #}

    if ($positives < 0) { $positives = '-'; }
    if ($gaps < 0) { $gaps = '-'; }

    if ($qchanged)
    {
        $qname =~ s/\s+$//;
        print "\n\@$qname\n$qlen\n";
        $qchanged = 0;
        $schanged = 1;
    }

    my $printed=0;
    if ($schanged)
    {
    	$printed=1;
        $sname =~ s/\s{2,}/ /g;
        $sname =~ s/\s+$//;
        $sname =~ s/^\s+//;
        $sname =~ s/^lcl\|//;
        print ">$sname\t$slen";
        $schanged = 0;
    }
    else
    {
	    $hchanged = 0;
	    $frame = '';
	    return;
    }

    if(!$printed)
    {
	    die "not printed $sname $qname ..";
    }
#else
#{
#	    if($PRINTALL)
#	    {
#		print ">$sname\t$slen";
#	    }
#	    else
#	    {
#		    return;
#	    }
#    }
    print "\t$ident\t$gaps\t$bits\t$expect\t$qstart\t$qend\t$sstart\t$send\t$HStrand\t$TStrand\n";
    $hchanged = 0;
    $frame = '';
    $TStrand=$HStrand="";

    $n_hits++;
}

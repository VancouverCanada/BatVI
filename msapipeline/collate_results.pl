#!/usr/bin/perl

use strict;
use warnings;

my $Lib=$ARGV[0];
open MSA,"$ARGV[1]" or die;
open PRED,"$ARGV[2]" or die;

#chr10   42385679        -       418     +       REPEAT  chr10:42385679:42385679:0
while(my $L=<MSA>)
{
	chomp $L;
	my ($Chr,$Loc,$SH,$Vir,$SV,$Type,$Contig)=split(/\t/,$L);
	print "$Lib\t$Type\t$Chr\t$Loc\t$Chr:$Loc:$SH\t$Vir$SV\t$Contig\n";
}

while(my $L=<PRED>)
{
	chomp $L;
	my ($SH,$Chr,$Loc,$L2,$Rank,$Rest)=split(/\t/,$L,6);
	if($Rank == 1)
	{
		print "$Lib\tWGS_RANK1\t$Chr\t$Loc\t$Chr:$Loc:$SH\n";
	}
}

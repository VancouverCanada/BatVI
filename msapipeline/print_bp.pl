#!/usr/bin/perl
use strict;
use warnings;

open BLAST,"$ARGV[0]" or die;

Loop:while(my $L=<BLAST>)
{
	if(substr($L,0,1) eq '@')
	{
		chomp $L;
		my $ID=substr($L,1);
		my $Contig_Len=<BLAST>;chomp $Contig_Len;
		my $Info=<BLAST>;chomp $Info;
		my($Genome,$Genome_length,$Ident,$Gaps,$Bits,$Expect,$Qstart,$Qend,$Sstart,$Send,$Plus,$Sign)=split(/\t/,$Info);
		$Genome=substr($Genome,1);
		die "Bad hit info ... $Info" if($Plus ne "+");
		
		my ($Multiplicity,$Chr,$Loc,$SignH,$Side,$Contig)=split(":",$ID,6);
		if($Side eq "L")
		{
			print "$Chr\t$Loc\t$SignH\t$Sstart\t$Sign\t$Multiplicity\t$Contig\n";
		}
		else
		{
			$SignH=&Flip_Sign($SignH);
			$Sign=&Flip_Sign($Sign);
			print "$Chr\t$Loc\t$SignH\t$Send\t$Sign\t$Multiplicity\t$Contig\n";

		}

	}
}

sub Flip_Sign()
{
	my $Sign=shift @_;
	if($Sign eq "-")
	{
		$Sign="+";
	}
	else
	{
		$Sign="-";
	}
	return $Sign;
}


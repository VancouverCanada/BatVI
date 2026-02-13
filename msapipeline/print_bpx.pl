#!/usr/bin/perl
use strict;
use warnings;

open BLAST,"$ARGV[0]" or die;
open FA,"$ARGV[1]" or die;

my %Sequences=();
Loop:while(my $L=<BLAST>)
{
	if(substr($L,0,1) eq '@')
	{
		chomp $L;
		my $V_Key=substr($L,1);
		my $Contig_Len=<BLAST>;chomp $Contig_Len;
		my $Info=<BLAST>;chomp $Info;
		my($Genome,$Genome_length,$Ident,$Gaps,$Bits,$Expect,$Qstart,$Qend,$Sstart,$Send,$Plus,$Sign)=split(/\t/,$Info);
		$Genome=substr($Genome,1);
		die "Bad hit info ... $Info" if($Plus ne "+");
		$Sequences{$V_Key}="$Sstart\t$Send\t$Sign";
		
	}
}

while(my $ID=<FA>)
{
	die "Bad fasta format.." if(substr($ID,0,1) ne '>');
	chomp $ID;
	$ID=substr($ID,1);
	my $V_Key=<FA>;die "Empty Virus Pattern..$ID" if(!$V_Key);chomp $V_Key;
	if(exists $Sequences{$V_Key})
	{
		my($Sstart,$Send,$Sign)=split(/\t/,$Sequences{$V_Key});


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
else
{
print STDERR "$ID-$V_Key\n"
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


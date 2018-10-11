#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  combine_hits.pl
#
#        USAGE:  ./combine_hits.pl LIB Max_Gap 
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Chandana Tennakoon (mn), mehner@fh-swf.de
#      COMPANY:  FH Südwestfalen, Iserlohn
#      VERSION:  1.0
#      CREATED:  10/29/2014 02:41:24 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

my %Chr_HashP=();
my %Chr_HashM=();
my $Merge_Gap=$ARGV[1];
my $NO_MSA=0;

my $FName=substr($ARGV[0],rindex($ARGV[0],"/")+1);
open PRED,"$ARGV[0]/predictions.opt.subopt.txt" or die;
open MSA,"$ARGV[0]/tmp.batvi/$FName.predictionsx.msa.txt" or $NO_MSA=1;

while(my $L=<PRED>)
{
	chomp $L;
	my @Fields=split(/\t/,$L);
	my $Sign=$Fields[0];
	my $Chr=$Fields[1];
	my $Loc=$Fields[2];
	my $VSign=$Fields[4];
	my $VLoc=$Fields[5];
	my $Count=$Fields[8];
	my $Rest=join("\t",@Fields[8..$#Fields]);
	if($Fields[7]==1)#Rank 1
	{
		if($Sign eq '+')
		{
			push(@{$Chr_HashP{$Chr}},join("\t",$Loc,$Sign,$VSign,$VLoc,$Count,$Rest));
		}
		else
		{
			die "Bad sign.." if($Sign ne '-');
			push(@{$Chr_HashM{$Chr}},join("\t",$Loc,$Sign,$VSign,$VLoc,$Count,$Rest));
		}
	}
}


my %Tophit=();
if(!$NO_MSA)
{
	while(my $L=<MSA>)
	{
		chomp $L;
		my @Fields=split(/\t/,$L);
		my $Sign=$Fields[2];
		my $Chr=$Fields[0];
		my $Loc=$Fields[1];
		my $VSign=$Fields[4];
		my $VLoc=$Fields[3];
		my $Type=$Fields[5];
		if($Type eq "TOPHIT")
		{
			if(!(exists $Tophit{"$Sign:$Chr:$Loc:$VSign:$VLoc"}))
			{
				$Tophit{"$Sign:$Chr:$Loc:$VSign:$VLoc"}=1;
			}	
		}
	}
}

my %Tophit_Hash=();

foreach my $K (keys %Tophit)
{
	my ($Sign,$Chr,$Loc,$VSign,$VLoc)=split(":",$K);
	push(@{$Tophit_Hash{$Chr}},join("\t",$Loc,$Sign,$VSign,$VLoc));
}

my @Final_Hits=();

foreach my $K (keys %Chr_HashP)
{
	my @A=@{$Chr_HashP{$K}};
	my($LocC,$SignC,$VSignC,$VLocC,$CountC,$ReadsC,$RestC)=split(/\t/,$A[0],7);
	for(my $i=1;$i<=$#A;$i++)
	{
		my($Loc,$Sign,$VSign,$VLoc,$Count,$Reads,$Rest)=split(/\t/,$A[$i],7);
		die "Unsorted Array..$A[0]" if($Loc<$LocC);
		die "Bad Sign.." if($SignC ne '+');

		if($LocC+$Merge_Gap>$Loc)#Nearby hits..
		{
			if($Reads>=$ReadsC)
			{
				($LocC,$SignC,$VSignC,$VLocC,$CountC,$ReadsC,$RestC)=($Loc,$Sign,$VSign,$VLoc,$Count,$Reads,$Rest);
			}
		}
		else
		{
			my $MSA_Hit=&Search_MSA($K,$LocC,$SignC);
			if($MSA_Hit)
			{
				push(@Final_Hits,"$MSA_Hit\t$ReadsC\t$RestC");
			}
			else
			{
				push(@Final_Hits,join("\t",$K,$LocC,$SignC,$VSignC,$VLocC,$CountC,$ReadsC,$RestC));
			}
			($LocC,$SignC,$VSignC,$VLocC,$CountC,$ReadsC,$RestC)=($Loc,$Sign,$VSign,$VLoc,$Count,$Reads,$Rest);
		}
	}
#Check Last/Only hit..
	my $MSA_Hit=&Search_MSA($K,$LocC,$SignC);
	if($MSA_Hit)
	{
		push(@Final_Hits,"$MSA_Hit\t$ReadsC\t$RestC");
	}
	else
	{
		push(@Final_Hits,join("\t",$K,$LocC,$SignC,$VSignC,$VLocC,$CountC,$ReadsC,$RestC));
	}

}

foreach my $K (keys %Chr_HashM)
{
	my @A=@{$Chr_HashM{$K}};
	my($LocC,$SignC,$VSignC,$VLocC,$CountC,$ReadsC,$RestC)=split(/\t/,$A[0],7);
	for(my $i=1;$i<=$#A;$i++)
	{
		my($Loc,$Sign,$VSign,$VLoc,$Count,$Reads,$Rest)=split(/\t/,$A[$i],7);
		die "Unsorted Array.." if($Loc<$LocC);
		die "Bad Sign.." if($SignC ne '-');

		if($LocC+$Merge_Gap>$Loc)#Nearby hits..
		{
			if($Reads>$ReadsC)
			{
				($LocC,$SignC,$VSignC,$VLocC,$CountC,$ReadsC,$RestC)=($Loc,$Sign,$VSign,$VLoc,$Count,$Reads,$Rest);
			}
		}
		else
		{
			my $MSA_Hit=&Search_MSA($K,$LocC,$SignC);
			if($MSA_Hit)
			{
				push(@Final_Hits,"$MSA_Hit\t$ReadsC\t$RestC");
			}
			else
			{
				push(@Final_Hits,join("\t",$K,$LocC,$SignC,$VSignC,$VLocC,$CountC,$ReadsC,$RestC));
			}
			($LocC,$SignC,$VSignC,$VLocC,$CountC,$ReadsC,$RestC)=($Loc,$Sign,$VSign,$VLoc,$Count,$Reads,$Rest);
		}
	}

	my $MSA_Hit=&Search_MSA($K,$LocC,$SignC);
	if($MSA_Hit)
	{
		push(@Final_Hits,"$MSA_Hit\t$ReadsC\t$RestC");
	}
	else
	{
		push(@Final_Hits,join("\t",$K,$LocC,$SignC,$VSignC,$VLocC,$CountC,$ReadsC,$RestC));
	}

}

print "LIB\tChr\tHuman Pos\tSign\tViral Sign\tViral Pos\tRead Count\tSplit Reads\tUniquely Mapped Reads\tMultiply Mapped Reads\tRank1 Hits\n"; 
foreach my $K (@Final_Hits)
{
	print "$ARGV[0]\t$K\n";
}

foreach my $Chr (keys %Tophit_Hash)
{
   my @A=@{$Tophit_Hash{$Chr}};
   L: for(my $i=0;$i<=$#A;$i++)
   {
	   if(!$A[$i]) {next L;}
	   my ($LocT,$SignT,$VSignT,$VLocT)=split(/\t/,$A[$i]);
	   print "$ARGV[0]\t".join("\t",$Chr,$LocT,$SignT,$VSignT,$VLocT,"MSA")."\n";
   }
}



sub Search_MSA()
{
	my $Chr=shift @_;my $Loc=shift @_;my $Sign=shift @_;
	my $MSA_Hit="";
	if($NO_MSA)
	{
		return $MSA_Hit;
	}

	if(exists $Tophit_Hash{$Chr})
	{
		my @A=@{$Tophit_Hash{$Chr}};
		L: for(my $i=0;$i<=$#A;$i++)
		{
			if(!$A[$i]) {next L;}
			my ($LocT,$SignT,$VSignT,$VLocT)=split(/\t/,$A[$i]);
			if($SignT eq $Sign)
			{
				if(abs($LocT-$Loc)<$Merge_Gap)
				{
					$MSA_Hit=join("\t",$Chr,$LocT,$SignT,$VSignT,$VLocT,"MSA");
					$A[$i]="";
					@{$Tophit_Hash{$Chr}}=@A;
					last L;
				}
			}
		}
	}
	return $MSA_Hit;
}

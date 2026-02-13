#!/usr/bin/perl
use strict;
use warnings;

open CONTIG,"$ARGV[0]" or die;
open BLAST,"$ARGV[1]" or die;

my %Contigs=();

my $L2;
while(my $L=<CONTIG>)
{
	die "Bad contig.." if(!($L2=<CONTIG>));
	chomp $L;chomp $L2;
	$L=substr($L,1);
	$Contigs{$L}=$L2;
}

Loop:while(my $L=<BLAST>)
{
	if(substr($L,0,1) eq '@')
	{
		chomp $L;
		my $ID=substr($L,1);
		my $Contig_Len=<BLAST>;chomp $Contig_Len;
		my $Contig=$Contigs{$ID};
		die "Incompatible contigs $Contig_Len - $Contig" if(length($Contig) != $Contig_Len);
		my $Info=<BLAST>;chomp $Info;
		my $Top_Expect=0;
		my @Hits=();
		my $Hit_Count=0;
		while($Info)
		{
			my($Genome,$Genome_length,$Ident,$Gaps,$Bits,$Expect,$Qstart,$Qend,$Sstart,$Send,$Plus,$Sign)=split(/\t/,$Info);
			$Top_Expect=$Expect if(!$Top_Expect);
			$Genome=substr($Genome,1);
			die "Bad hit info ... $Info" if($Plus ne "+");
			if($Expect eq $Top_Expect)
			{
				$Hit_Count++;
				if($Qstart>=1)# case 1) .......|HG
				{
					if($Qstart>5)#avoid very small virus..
					{
						my $Vir_Clip;
						if($Qstart>25)
						{
							my $T=$Qstart-25;
							$Vir_Clip=substr($Contig,$Qstart-25,25);
						}
						else
						{
							$Vir_Clip=substr($Contig,0,$Qstart);
						}
						if($Vir_Clip)
						{
							push(@Hits,"$Genome:$Sstart:$Sign:R:$ID\n$Vir_Clip\n");
						}
						else
						{
							print $Info;
						}
					}
				}
				if($Qend<$Contig_Len)	# case 2) HG|....
				{
					my $Vir_Clip;
					if(($Contig_Len-$Qend)>25)
					{
						$Vir_Clip=substr($Contig,$Qend,25);
					}
					else
					{
						$Vir_Clip=substr($Contig,$Qend);
					}
					push (@Hits,"$Genome:$Send:$Sign:L:$ID\n$Vir_Clip\n");
				}
			}
			$Info=<BLAST>;chomp $Info if($Info);
		}
		foreach my $Line (@Hits)
		{
			if($Hit_Count == 1)#Uniq top hit..
			{
				print ">TOPHIT:$Line"; 
			}
			else	
			{
				print ">REPEAT:$Line"; 
			}
		}


	}
}

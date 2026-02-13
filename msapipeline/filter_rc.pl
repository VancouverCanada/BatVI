#!/usr/bin/perl
#FIlters the doubly entered msa RC/+ pairs
use warnings;
use strict;

my $Last_ID="";
my @Reads_Array=();
my $Count;
while (my $L=<>)
{
	chomp $L;
	my $Read=<>;
	@_=split(":",$L);
	$L="$_[0]:$_[1]:$_[2]";
	if($Last_ID ne $L)
	{
		for(my $i=0;$i<=$#Reads_Array;$i++)
		{
			if($Reads_Array[$i])
			{
				my $RC=&Reverse_Complement($Reads_Array[$i]);
				Check_Loop:for(my $j=0;$j<=$#Reads_Array;$j++)
				{
					next Check_Loop if($i==$j);
					if($RC eq $Reads_Array[$j])
					{
						$Reads_Array[$j]=0;
						print "$Last_ID:$Count\n$RC\n";
						$Count++;
						last Check_Loop;
					}
				}
			}
		}
		$Last_ID=$L;
		@Reads_Array=();$Count=0;
	}
	chomp $Read;
	push(@Reads_Array,$Read);
}

for(my $i=0;$i<=$#Reads_Array;$i++)
{
	if($Reads_Array[$i])
	{
		my $RC=&Reverse_Complement($Reads_Array[$i]);
Check_Loop2:for(my $j=0;$j<=$#Reads_Array;$j++)
	   {
		   next Check_Loop2 if($i==$j);
		   if($RC eq $Reads_Array[$j])
		   {
			   $Reads_Array[$j]=0;
			   print "$Last_ID:$Count\n$RC\n";
			   $Count++;
			   last Check_Loop2;
		   }
	   }
	}
}

sub Reverse_Complement()
{
	my $DNA = shift;
	my $RC = reverse($DNA);
	$RC =~ tr/ACGTacgt/TGCAtgca/;
	return $RC;
}

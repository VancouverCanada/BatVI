#!/usr/bin/perl
#Opens unmapped fasta read file (arg 1, arg2), find reads given in arg 1
#then looks for stdin for matching reads..
use strict;
use warnings;

open L,"$ARGV[0]/tmp/hbv.final.list" or die;
open U1,"$ARGV[0]/fastq/unmapped_1.fa" or die;
open U2,"$ARGV[0]/fastq/unmapped_2.fa" or die;
open HEAD,">$ARGV[0]/tmp/head.fq" or die;
open TAIL,">$ARGV[0]/tmp/tail.fq" or die;

my %Des=();
my $Rev_Count=0;
while(my $L=<L>)
{
	chomp $L;
	$Des{$L}= "";
}

my $U2_Blank=0;
Loop:while(my $L=<U1>)
{
	my $L1=<U1>;
	if(substr($L,-3,1) eq "/")
	{
		$L=substr($L,1,(length(chomp $L)-3));
	}
	else
	{
		chomp $L;
		$L=substr($L,1);
	}

	if(exists $Des{$L})
	{
		if(!$Des{$L}) #Blank read..
		{
			$Des{$L}=$L1;
		}
		else
		{
			print HEAD ">$L\n$L1";
			print TAIL ">$L\n$Des{$L}";
			delete $Des{$L};
		}
	}

	if($U2_Blank)
	{
		next Loop;
	}
	else
	{
		if(my $L=<U2>)
		{
			my $L1=<U2>;
			if(substr($L,-3,1) eq "/")
			{
				$L=substr($L,1,(length(chomp $L)-3));
			}
			else
			{
				chomp $L;
				$L=substr($L,1);
			}
			if(exists $Des{$L})
			{
				if(!$Des{$L}) #Blank read..
				{
					$Des{$L}=$L1;
				}
				else
				{
					print TAIL ">$L\n$L1";
					print HEAD ">$L\n$Des{$L}";
					delete $Des{$L};
				}
			}
		}
		else
		{
			$U2_Blank=1;
		}
	}
	
}

if(!$U2_Blank)
{
	while (my $L=<U2>)
	{
		my $L1=<U2>;
		if(substr($L,-3,1) eq "/")
		{
			$L=substr($L,1,(length(chomp $L)-3));
		}
		else
		{
			chomp $L;
			$L=substr($L,1);
		}
		if(exists $Des{$L})
		{
			if(!$Des{$L}) #Blank read..
			{
				$Des{$L}=$L1;
			}
			else
			{
				print TAIL ">$L\n$L1";
				print HEAD ">$L\n$Des{$L}";
				delete $Des{$L};
			}
		}
	}
}

foreach my $Key (keys %Des)
{
	print "$Key Not found\n";
}
print STDERR "Reveresed reads: $Rev_Count\n";
exit;

while(my $L=<STDIN>)
{
	my @A=split(/\t/,$L);
	if(exists $Des{$A[0]})
	{
		if(!$Des{$A[0]})
		{
			die "$A[0] Not matched $Des{$A[0]}..";
		}
		else
		{
			if($A[1] & 16) #Reveresd read..
			{
				$A[9]=&Reverse_Complement($A[9]);	
				$Rev_Count++;
			}
			if($A[1] & 64)
			{
				print HEAD ">$A[0]\n$A[9]\n";
				print TAIL ">$A[0]\n$Des{$A[0]}";
			}
			elsif ($A[1] & 128)
			{
				print TAIL ">$A[0]\n$A[9]\n";
				print HEAD ">$A[0]\n$Des{$A[0]}";
			}
			else
			{
				die "$A[0] Flag error matched..";
			}
			delete $Des{$A[0]};
		}
		if(!%Des)#All keys processed..
		{
			exit;
		}
	}
}

foreach my $Key (keys %Des)
{
	print "$Key Not found\n";
}
print STDERR "Reveresed reads: $Rev_Count\n";

sub Reverse_Complement()
{
	my $DNA = shift;
	my $RC = reverse($DNA);
	$RC =~ tr/ACGTacgt/TGCAtgca/;
	return $RC;
}

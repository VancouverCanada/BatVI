#!/usr/bin/perl
use strict;
use warnings;

my @Hit_Array;
my $Current_Read;
while(my $L=<>)
{
	if(substr($L,0,1) eq "@")
	{
		foreach my $K (@Hit_Array)
		{
			print "$K\t$Current_Read\n";
		}
		@Hit_Array=();
		chomp $L;
		$Current_Read=$L;
	}
	if(substr($L,0,1) eq '>')
	{
		@_=split(/\t/,$L);
		my $Loc=$_[8];my $Chr=substr($_[0],1);
		push(@Hit_Array,"$Chr\t$Loc");
	}
}

foreach my $K (@Hit_Array)
{
	print "$K\t$Current_Read\n";
}

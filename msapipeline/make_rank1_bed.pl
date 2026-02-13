#!/usr/bin/perl
use warnings;
use strict;

my $ChrL="";
my $LocL;
while(my $L=<>)
{
	my ($Chr,$Loc)=split(/\t/,$L);
	if($ChrL)
	{
		if($Chr eq $ChrL)
		{
			die "Unsorted file.." if($Loc<$LocL);
			if(($Loc-$LocL)<1000)
			{
				$LocL=$Loc;
			}
			else
			{
				my $St=$LocL-100;my $Ed=$LocL+100;
				print "$ChrL\t$St\t$Ed\n";
			}
		}
		else
		{
			my $St=$LocL-100;my $Ed=$LocL+100;
			print "$ChrL\t$St\t$Ed\n";
		}
		$ChrL=$Chr;$LocL=$Loc;
	}
	else
	{
			$ChrL=$Chr;$LocL=$Loc;
	}
}
my $St=$LocL-100;my $Ed=$LocL+100;
print "$ChrL\t$St\t$Ed\n";

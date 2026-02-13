#!/usr/bin/perl

$L1=<>;chomp $L1;
while($L2=<>)
{
	chomp $L2;
	my ($Chr1,$Loc1,$Sign1,$VLoc2,$VS1)=split(/\t/,$L1);
	my ($Chr2,$Loc2,$Sign2,$VLoc2,$VS2)=split(/\t/,$L2);
	if(($Chr1 eq $Chr2)&&($Sign1 ne $Sign2)&&(abs($Loc1-$Loc2)<1000))
	{
		if($VS1 eq $VS2)
		{
			print "$L1\tPAIRED\n";
			print "$L2\tPAIRED\n";
		}
	}
	$L1=$L2;
}

#!/usr/bin/perl
use strict;
use warnings;

my $L1;my $L2;
my $Parsed=0;
Loop:while($L1=<>)
{
	$Parsed=1;
	next Loop if(substr($L1,0,1) eq '@');
	my @Sam1=split(/\t/,$L1);
	die "Unpaired sam .. $L1" if(!($L2=<>));
	my @Sam2=split(/\t/,$L2);
	die "Unpaired SAM $Sam1[0]\t$Sam2[0]" if($Sam1[0] ne $Sam2[0]);
	my $Len1=length $Sam1[9];
	my $Len2=length $Sam2[9];
	die "Unpaired SAM length $Len1\t$Len2" if($Len1 ne $Len2);

	my $Match_Len1=0;
	my $Match_Len2=0;
	if($Sam1[5]=~/[a-zA-Z]*([0-9]*)M/)
	{
		$Match_Len1=$1;
	}
	if($Sam2[5]=~/[a-zA-Z]*([0-9]*)M/)
	{
		$Match_Len2=$1;
	}
	if(($Match_Len1 != $Len1)&&($Match_Len2 != $Len2))
	{
		print ">$Sam1[0]\n$Sam1[9]\n";
	}

}
if(!$Parsed)
{
	die "Empty stream..";
}

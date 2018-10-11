#!/usr/bin/perl

open FILE,$ARGV[0] or die;
open OUT,">$ARGV[1]" or die;
$Mark=$ARGV[2];#ID to show library name..
$Insert=$ARGV[3];
#$Key="$Mark:$Insert:";
$Key="";

$L1=<FILE>;
split(/\t/,$L1);
$Last=$_[0];

Loop:while($L2=<FILE>)
{
	split(/\t/,$L2);
	$Current=$_[0];
	if($Current eq $Last)#Match..
	{
		print OUT "$Key:$L1$Key:$L2";
		$L1=<FILE>;
		if(!$L1)
		{
			last Loop;
		}
		split(/\t/,$L1);
		$Last=$_[0];
	}
	else
	{
		print STDERR "$Key:$L1";
		$L1=$L2;
		$Last=$Current;
	}
}

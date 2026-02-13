#!/usr/bin/perl

while($L=<>)
{
	chomp $L;
	($C,$L1,$L2)=split(":",$L,3);
	print "$C\t$L1\t$L2\n";
}


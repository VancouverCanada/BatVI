#!/usr/bin/perl

$Real_Count=0;
while($L=<>)
{
	if(substr($L,0,1) eq '*')
	{
		die "Bad format.. $L" if(!($ID=<>));
		die "Bad format.. $L" if(!($L=<>));
		die "Bad format.. $L" if (substr($L,0,1) ne '*');
		chomp $ID;
		$ID=~ s/\..\.fa//;
		$Count=1;
	}
	else
	{
		die "Bad format.. $L" if (substr($L,0,1) ne '>');
		die "Bad format $L" if(!($L=<>));
		print ">$ID:$Count:$Real_Count\n$L";$Count++;$Real_Count++;
	}
}

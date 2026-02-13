#!/usr/bin/perl
$Count=0;
$Pass=0;
while($L=<>)
{
	chomp $L;
	$L= substr($L,rindex($L,"/"));#remove path..
	$L= substr($L,1,length($L)-6);#remove extension..
	$Pass=0 if($Pass==2);
	$Pass++;

	if($Pass==2)
	{
		die "Bad Format.." if($First_ID ne $L);
	}
	else
	{
		$First_ID=$L;
		print "$L\n";
	}
	$L=<>;chomp $L;

	while($L ne "--")
	{
		if(substr($L,0,1) eq '>')
		{
			print "$L:$Pass\n";
		}
		else
		{
			print "$L\n";
		}
		$L=<>;chomp $L;
	}
	print "$L\n" if($Pass==2);
}

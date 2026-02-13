#!/usr/bin/perl


open F,$ARGV[0] or die;
if(!$ARGV[1])
{
	$Clusterdist=1000;
}
else
{
	$Clusterdist=$ARGV[1];
}

$Cluster_Size=$ARGV[2];
$Last_Chr="";
$Last_Loc=0;
while($L=<F>)
{
	@_=split (/\t/,$L);
	$Chr=$_[0];$Loc=$_[1];
	if($Last_Chr ne $Chr)
	{
		$ASize=$#Cluster+1;
		if($Last_Chr && $ASize >=$Cluster_Size )#Not the first pass..
		{
			$ID="$Last_Chr:$St:$Ed";
			print STDERR "$ID:$ASize\n";
			foreach $Line (@Cluster)
			{
				print "$ID\t$Line";
			}
		}
		$Last_Chr=$Chr;
		$St=$Loc;$Ed=$Loc;
		$Last_Loc=$Loc;
		@Cluster=();push (@Cluster,$L);
	}	
	else
	{
		if($Last_Loc>$Loc)
		{
			die "$Last_Loc $Loc Unsorted file..";
		}
		elsif($Loc-$Last_Loc < $Clusterdist) #clustering distance...
		{
			if($Loc<$St) { $St=$Loc; }
			if($Loc>$Ed) { $Ed=$Loc; }
			push (@Cluster,$L);
		}
		else
		{
			if($Last_Chr)#Not the first pass..
			{	
				$ASize=$#Cluster+1;
				if($ASize >=$Cluster_Size)
				{
					$ID="$Last_Chr:$St:$Ed";
					print STDERR "$ID:$ASize\n";
					foreach $Line (@Cluster)
					{
						print "$ID\t$Line";
					}
				}
			}
			else
			{ 
				die "format error.."; 
			}
			$Last_Chr=$Chr;
			$St=$Loc;$Ed=$Loc;
			$Last_Loc=$Loc;
			@Cluster=();push (@Cluster,$L);
		}
		$Last_Loc=$Loc;
	}
}

if($Last_Chr)
{
	$ASize=$#Cluster+1;
	if($ASize >=$Cluster_Size)#Not the first pass..
	{
		$ID="$Last_Chr:$St:$Ed";
		print STDERR "$ID:$ASize\n";
		foreach $Line (@Cluster)
		{
			print "$ID\t$Line";
		}
	}
}
else
{ 
	die "format error.."; 
}

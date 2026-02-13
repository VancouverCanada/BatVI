#!/usr/bin/perl
use strict;
use warnings;
#chr1    614873  -       2365    +       REPEAT  chrUn_gl000240:28478:29112:0
my %Plus_Integrations=();
my %Plus_Tophit=();
my %Plus_IntegrationsV=();
my %Minus_Integrations=();
my %Minus_IntegrationsV=();
my %Minus_Tophit=();
my $FUSEGAP=10;
open FILE,"$ARGV[0]" or die;
while(my $L=<FILE>)#split integrations into two halves..
{
	chomp $L;
	my ($Chr,$LocH,$SignH,$LocV,$SignV,$Type,$Cluster)=split(/\t/,$L);
	if($SignH eq "+")
	{
		push(@{$Plus_Integrations{$Chr}},$LocH);
		push(@{$Plus_IntegrationsV{"$Chr:$LocH"}},"$SignV\t$LocV");
		if($Type eq "TOPHIT")
		{
			push(@{$Plus_Tophit{$Chr}},$LocH);
		}
	}
	else
	{
		push(@{$Minus_Integrations{$Chr}},$LocH);
		push(@{$Minus_IntegrationsV{"$Chr:$LocH"}},"$SignV\t$LocV");
		if($Type eq "TOPHIT")
		{
			push(@{$Minus_Tophit{$Chr}},$LocH);
		}
	}
}

#########################################
#Collapse Plus Integrations
#########################################

my %Plus_Integrations_Collapsed=();
my %Virus_Side=();
foreach my $Chr (keys %Plus_Integrations)
{
	my @Hits= sort { $a <=> $b } @{$Plus_Integrations{$Chr}};
	my %Collapsible_Hits=();
	push(@{$Collapsible_Hits{$Hits[0]}},$Hits[0]);
	my $Last_Hit=$Hits[0];
	my %Vir_Side=();#Store viral contig info..

	
	for (my $i=1;$i<=$#Hits;$i++)
	{
		my $Dist=($Hits[$i]-$Last_Hit);
		die "Bad co-ordinated $Hits[$i]\t$Last_Hit" if($Dist<0);
		if($Dist<$FUSEGAP)
		{
			push(@{$Collapsible_Hits{$Last_Hit}},$Hits[$i]);
			push(@{$Vir_Side{"$Chr:$Last_Hit"}},@{$Plus_IntegrationsV{"$Chr:$Hits[$i]"}});
		}
		else
		{
			push(@{$Vir_Side{"$Chr:$Hits[$i]"}},@{$Plus_IntegrationsV{"$Chr:$Hits[$i]"}});
			push(@{$Collapsible_Hits{$Hits[$i]}},$Hits[$i]);
		}
	}

	foreach my $Collapsed_Hits (keys %Collapsible_Hits)
	{
		my @Hit_List=@{$Collapsible_Hits{$Collapsed_Hits}};
		my $Rightmost=$Hit_List[0];
		for(my $i=1;$i<=$#Hit_List;$i++)
		{
			if($Hit_List[$i]>$Rightmost)
			{
				$Rightmost=$Hit_List[$i];
			}
		}
		push(@{$Plus_Integrations_Collapsed{$Chr}},$Rightmost);
		
		@Hit_List=@{$Vir_Side{"$Chr:$Collapsed_Hits"}};

		my $S;my $Old_S="";my $L;my $Leftmost; 
		foreach my $VHits (@Hit_List)
		{
			($S,$L)=split(/\t/,$VHits);
			if(!$Old_S){$Old_S=$S;$Leftmost=$L;};
			if($S ne $Old_S){$L="?",$S="?"}
			else
			{
				if($L<$Leftmost)
				{
					$Leftmost=$L;
				}
			}
		}	
		$Virus_Side{"$Chr:$Rightmost:+"}="$Leftmost:$S";
	}
}

#########################################
#Collapse Minus Integrations
#########################################

my %Minus_Integrations_Collapsed=();
foreach my $Chr (keys %Minus_Integrations)
{
	my @Hits= sort { $a <=> $b } @{$Minus_Integrations{$Chr}};
	my %Collapsible_Hits=();
	push(@{$Collapsible_Hits{$Hits[0]}},$Hits[0]);
	my $Last_Hit=$Hits[0];
	my %Vir_Side=();#Store viral contig info..

	
	for (my $i=1;$i<=$#Hits;$i++)
	{
		my $Dist=($Hits[$i]-$Last_Hit);
		die "Bad co-ordinated $Hits[$i]\t$Last_Hit" if($Dist<0);
		if($Dist<$FUSEGAP)
		{
			push(@{$Collapsible_Hits{$Last_Hit}},$Hits[$i]);
			push(@{$Vir_Side{"$Chr:$Last_Hit"}},@{$Minus_IntegrationsV{"$Chr:$Hits[$i]"}});
		}
		else
		{
			push(@{$Vir_Side{"$Chr:$Hits[$i]"}},@{$Minus_IntegrationsV{"$Chr:$Hits[$i]"}});
			push(@{$Collapsible_Hits{$Hits[$i]}},$Hits[$i]);
		}
	}

	foreach my $Collapsed_Hits (keys %Collapsible_Hits)
	{
		my @Hit_List=@{$Collapsible_Hits{$Collapsed_Hits}};
		my $Leftmost=$Hit_List[0];
		for(my $i=1;$i<=$#Hit_List;$i++)
		{
			if($Hit_List[$i]<$Leftmost)
			{
				$Leftmost=$Hit_List[$i];
			}
		}
		push(@{$Minus_Integrations_Collapsed{$Chr}},$Leftmost);

		@Hit_List=@{$Vir_Side{"$Chr:$Collapsed_Hits"}};

		my $S;my $Old_S="";my $L;my $Rightmost; 
		foreach my $VHits (@Hit_List)
		{
			($S,$L)=split(/\t/,$VHits);
			if(!$Old_S){$Old_S=$S;$Rightmost=$L;};
			if($S ne $Old_S){$L="?",$S="?"}
			else
			{
				if($L<$Rightmost)
				{
					$Rightmost=$L;
				}
			}
		}	
		$Virus_Side{"$Chr:$Leftmost:-"}="$Rightmost:$S";
	}
}
##############################################
# convert tophits to collapsed hits
##############################################

my %TophitU=();#Keeps the final list of tophits
my %Tophit_ConflictU=();#Keeps the final list of conflicting tophits
my %RepeatU=();#Keeps the final list of repeats 
my %Repeat_FalseU=();#Keeps the final list of false repeats 

my %Plus_TophitU=();
foreach my $Chr (keys %Plus_Tophit)
{
	my @A=@{$Plus_Tophit{$Chr}};
	for(my $i=0;$i<=$#A;$i++)
	{
		my @P=@{$Plus_Integrations_Collapsed{$Chr}};
		my $Nearest=&Locate_Nearest($A[$i],\@P,20);
		die "Missing intrgration $Chr:$A[$i]" if($Nearest== -1);
		$A[$i]=$Nearest;
	}
	my %Tmp=();
	$Plus_Tophit{$Chr} = [ grep { ! $Tmp{$_} ++ } @A ];
	foreach my $Loc (@{$Plus_Tophit{$Chr}})
	{
		$Plus_TophitU{"$Chr:$Loc"}=1;
		$TophitU{"$Chr:$Loc:+"}=1;
	}
}



my %Minus_TophitU=();
foreach my $Chr (keys %Minus_Tophit)
{
	my @A=@{$Minus_Tophit{$Chr}};
	for(my $i=0;$i<=$#A;$i++)
	{
		my @P=@{$Minus_Integrations_Collapsed{$Chr}};
		my $Nearest=&Locate_Nearest($A[$i],\@P,20);
		die "Missing intrgration $Chr:$A[$i]" if($Nearest== -1);
		$A[$i]=$Nearest;
	}
	my %Tmp=();
	$Minus_Tophit{$Chr} = [ grep { ! $Tmp{$_} ++ } @A ];
	foreach my $Loc (@{$Minus_Tophit{$Chr}})
	{
		$Minus_TophitU{"$Chr:$Loc"}=1;
		$TophitU{"$Chr:$Loc:-"}=1;
	}
}

seek FILE, 0, 0;
my $Last_Cluster="";
open REPLOG,">$ARGV[1]" or die;
my %Conflicting_Tophits=();

my @Tophit_Array=();
my @Repeat_Array=();
while(my $L=<FILE>)
{
	chomp $L;
	my ($Chr,$LocH,$SignH,$LocV,$SignV,$Type,$Cluster)=split(/\t/,$L);
	if($Last_Cluster ne $Cluster)
	{
		if($Last_Cluster)#Ignore first pass..
		{
			my $Tophit_Count=scalar(@Tophit_Array);
			print REPLOG "\@$Last_Cluster\t$Tophit_Count\n";
			my $Tophit_Present=0;
			if($Tophit_Count == 1)#Only one tophit
			{
				$Tophit_Present=1;
				print REPLOG "$Tophit_Array[0]\tTOPHIT\n";
			}
			elsif($Tophit_Count >1)#Many tophits
			{
				my @Tophits_Resolved=&Resolve_Tophit(\@Tophit_Array);
				if(@Tophits_Resolved)
				{
					$Tophit_Present=1;
					for(my $i=0;$i<=$#Tophits_Resolved;$i++)
					{
						print REPLOG "$Tophits_Resolved[$i]\tTOPHIT\n";
					}
				}
				else
				{
					for(my $i=0;$i<=$#Tophit_Array;$i++)
					{
						delete $TophitU{$Tophit_Array[$i]}if(exists $TophitU{$Tophit_Array[$i]});
						$Tophit_ConflictU{$Tophit_Array[$i]}=1;
						print REPLOG "$Tophit_Array[$i]\tTOPHIT_CONFLICT\n";
					}
				}
			}
			foreach my $K (@Repeat_Array)
			{
				if($Tophit_Present)
				{
					$Repeat_FalseU{$K}=1;
					print REPLOG "$K\tFALSE_REPEAT\n";
				}
				else
				{
					$RepeatU{$K}=1;
					print REPLOG "$K\tREPEAT\n";
				}
			}
		}
		@Tophit_Array=();
		@Repeat_Array=();
		$Last_Cluster=$Cluster;
	}

	my $Tophit=0;
	if($SignH eq "+")
	{
#Substitute unified BP..
		my @P=@{$Plus_Integrations_Collapsed{$Chr}};
		my $LocH=&Locate_Nearest($LocH,\@P,20);
		die "Missing intrgration $Chr:$LocH" if($LocH== -1);
		$Tophit=1 if(exists $Plus_TophitU{"$Chr:$LocH"});

	}
	else
	{
		die "Bad sign.." if($SignH ne "-");
#Substitute unified BP..
		my @P=@{$Minus_Integrations_Collapsed{$Chr}};
		my $LocH=&Locate_Nearest($LocH,\@P,20);
		die "Missing intrgration $Chr:$LocH" if($LocH== -1);
		$Tophit=1 if(exists $Minus_TophitU{"$Chr:$LocH"});
	}
	if($Tophit)
	{
		push(@Tophit_Array,"$Chr:$LocH:$SignH");
	}
	else
	{
		push(@Repeat_Array,"$Chr:$LocH:$SignH");
	}
}

###############################
#Print TOPHITS
###############################
#foreach my $K (keys %Plus_TophitU)
#{
#	print "$K\tTOPHITP\n";
#}
#foreach my $K (keys %Minus_TophitU)
#{
#	print "$K\tTOPHITM\n";
#}

foreach my $K (keys %TophitU)
{
	@_=split(":",$K);my $K1="$_[0]\t$_[1]\t$_[2]";
	@_=split(":",$Virus_Side{$K});my $K2="$_[0]\t$_[1]";
	print "$K1\t$K2\tTOPHIT\n";
}

foreach my $K (keys %Tophit_ConflictU)
{
	@_=split(":",$K);my $K1="$_[0]\t$_[1]\t$_[2]";
	@_=split(":",$Virus_Side{$K});my $K2="$_[0]\t$_[1]";
	print "$K1\t$K2\tTOPHAT_CONFLICT\n";
}

foreach my $K (keys %Repeat_FalseU)
{
	delete $RepeatU{$K} if(exists $RepeatU{$K});
	@_=split(":",$K);my $K1="$_[0]\t$_[1]\t$_[2]";
	@_=split(":",$Virus_Side{$K});my $K2="$_[0]\t$_[1]";
	print "$K1\t$K2\tREPEAT_FALSE\n";
}

foreach my $K (keys %RepeatU)
{
	@_=split(":",$K);my $K1="$_[0]\t$_[1]\t$_[2]";
	@_=split(":",$Virus_Side{$K});my $K2="$_[0]\t$_[1]";
	print "$K1\t$K2\tREPEAT\n";
}

sub Resolve_Tophit()
{
	my @A = @{$_[0]};
	if(scalar(@A)!=2)
	{
		@A=();
	}
	else
	{
		my ($Chr1,$Loc1,$Sign1)=split(":",$A[0]);
		my ($Chr2,$Loc2,$Sign2)=split(":",$A[1]);
		if(($Chr1 ne $Chr2)||($Sign1 eq $Sign2)||(abs($Loc1-$Loc2)>1000))
		{
			@A=();
		}
	}
	return @A;
}

sub Locate_Nearest()
{
	my $Value=$_[0];
	my @A = @{$_[1]};
	my $Dist=$_[2];
	Locate:foreach my $K (@A)
	{
		if(abs($K-$Value)<=$Dist)
		{
			return $K; 
		}
	}
       return -1;
}

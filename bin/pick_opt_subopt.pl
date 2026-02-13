#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  pick_opt_subopt.pl
#
#        USAGE:  ./pick_opt_subopt.pl Blast_File 
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Chandana Tennakoon (mn), mehner@fh-swf.de
#      CREATED:  10/04/2014 04:49:32 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

my @Top_Hit=();my @Sub_Opt=();
Loop: while(my $L=<>)
{
	if(substr($L,0,1) eq '@')
	{
		@Top_Hit=();@Sub_Opt=();
		print $L;$L=<>;print $L;
#$L=<>;
#@_=split(/\t/,$L,7);
		my $Top=-1;my $SubOpt= 2;
		my $Pass=0;
Loop2:  while(!$Pass)
		{
			$L=<>;
			last Loop2 if eof;
			if(substr($L,0,1) eq '>') 
			{
				@_=split(/\t/,$L,7);
				my $Exp=$_[5];
				$Top=$Exp if($Top==-1);
				if($Exp eq $Top)
				{
					push (@Top_Hit,$L);
				}
				elsif ($Exp eq $SubOpt)
				{
					push (@Sub_Opt,$L);
				}
				elsif($Exp<$SubOpt)
				{
					@Sub_Opt=();
					push (@Sub_Opt,$L);
					$SubOpt=$Exp;
				}
			}
			else #Blank..
			{
				$Pass=1;#Break loop..
			}
		}
	}
	foreach my $Line (@Top_Hit)
	{
		print $Line;
	}
if(@Top_Hit==1)
	{
		foreach my $Line (@Sub_Opt)
		{
			print $Line;
		}
	}
	print $L if($L eq "\n");
	last Loop if eof;
}



#!/usr/bin/perl

use warnings;
use strict;

open LIST,"$ARGV[0]" or die;
open H,"$ARGV[1]" or die;
open T,"$ARGV[2]" or die;

my %File_Handle_H=();
my %File_Handle_T=();
my %Reads=();


while(my $L=<LIST>)
{
	chomp $L;
	my ($File,$ID)=split(/\t/,$L);
	$ID=substr($ID,1);
	if(!exists $File_Handle_H{$File})
	{
		open $File_Handle_H{$File},">reads/$File.1.fa" or die;
		open $File_Handle_T{$File},">reads/$File.2.fa" or die;
	} 
	$Reads{$ID}=$File;
}

my $H1;my $H2;my $T1;my $T2;
while($H1=<H>)
{
	die "Read mismatch" if (!($H2=<H>));
	die "Read mismatch" if (!($T1=<T>));
	die "Read mismatch" if (!($T2=<T>));

	chomp $H1;chomp $H2;
	chomp $T1;chomp $T2;
	my $ID=substr($H1,1);
	if(exists $Reads{$ID})
	{
		print {$File_Handle_H{$Reads{$ID}}} "$H1\n$H2\n";
		print {$File_Handle_T{$Reads{$ID}}} "$T1\n$T2\n";
	}
}

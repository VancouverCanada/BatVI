#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  split.pl
#
#        USAGE:  ./split.pl in_file out_file_pref number_of_files  
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Chandana Tennakoon (mn), mehner@fh-swf.de
#      COMPANY:  FH Südwestfalen, Iserlohn
#      VERSION:  1.0
#      CREATED:  10/08/2014 02:52:17 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

my $N=$ARGV[2];
my $Pref=$ARGV[1];
open IN,"$ARGV[0]" or die("$!");

my @Split_Files=();
for(my $i=0;$i<$N;$i++)
{
	open $Split_Files[$i],">${Pref}.$i" or die("$!");
}

my $Count= -1;
while(my $L=<IN>)
{
	$Count++;
	$Count=0 if($Count==$N);
	print {$Split_Files[$Count]} $L;
}



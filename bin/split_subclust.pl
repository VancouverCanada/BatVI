#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  split_subclust.pl
#
#        USAGE:  ./split_subclust.pl  
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
#      CREATED:  11/09/2014 08:06:14 AM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

while (my $L=<>)
{
	@_=split(/\t/,$L);
	my ($Chr,$St,$Ed)=split(":",$_[0]);
	while($St<$Ed)
	{
		print "$Chr\t$St\n";$St+=270;
	}
	print "$Chr\t$Ed\n";
}


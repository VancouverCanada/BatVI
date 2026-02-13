#!/usr/bin/perl
use strict;
use warnings;

my $Count=0;


my $directory = $ARGV[0];
opendir (DIR, $directory) or die $!;

for my $file (sort readdir(DIR)) 
{

# We only want *.fa
  next unless (-f "$ARGV[0]/$file");
  next unless ($file =~ m/\.fa$/);
        
  print STDERR "Processing $file\n";
  open F,"$ARGV[0]/$file" or die $!;
  
  print "$ARGV[0]/$file\n";
  while(my $L=<F>)
  {
    print $L;
    if(substr($L,0,1) ne '>')
    {
      chomp $L;
      my $RC=&Reverse_Complement($L);
      print ">$Count\n$RC\n";
      $Count++;
    }
  }
  print "--\n";
  close F;
}

sub Reverse_Complement()
{
	my $DNA = shift;
	my $RC = reverse($DNA);
	$RC =~ tr/ACGTacgt/TGCAtgca/;
	return $RC;
}


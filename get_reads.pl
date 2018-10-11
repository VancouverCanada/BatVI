#!/usr/bin/perl

use strict;
use warnings;

if (!@ARGV)
{
  print STDERR "get_reads.pl directory Chr Start [-c|-x]\n";
  print STDERR "Picks possible reads containing integrations at Chr:Start in the directory\n";
  print STDERR "-c will optionally output the cluster ID used by batvi\n";
  print STDERR "-x will extract the actual reads, otherwise only the read ID will be printed\n";
  die;
}

open DATA, "$ARGV[0]/t.opt.subopt.cluster" or die;
my $Chr=$ARGV[1];
my $St=$ARGV[2];
my $Cluster=0;
my $Extract=0;
if (exists $ARGV[3])
{
  $Cluster=1 if($ARGV[3] eq "-c");
  $Extract=1 if($ARGV[3] eq "-x");
}
my %Hash=();

while(my $L=<DATA>)
{
  chomp $L;
  my @F=split(' ',$L);
  my @L=split(":",$F[0]);
  if($L[0] eq $Chr)
  {
    if(($L[1] <= $St) && ($L[2] >= $St))
    {
      my $ID=substr($F[3],1);
      $ID =~ s/>//;
      $Hash{$ID}=1;
      if(!$Extract)
      {
        print $ID;
        print $F[0] if($Cluster);
        print "\n";
      }
    }
  }
}

if($Extract)
{
  open H, "$ARGV[0]/fastq/head.fa" or die;
  open T, "$ARGV[0]/fastq/tail.fa" or die;
  while(my $H=<H>)
  {
    my $H_seq=<H>;
    my $T=<T>;
    my $T_seq=<T>;

    chomp $H;
    my $ID=substr($H,1);
    $ID =~ s/>//;
    if(exists $Hash{$ID})
    {
      print ">${ID}_1\n$H_seq";
      print ">${ID}_2\n$T_seq";
    }
  }
}

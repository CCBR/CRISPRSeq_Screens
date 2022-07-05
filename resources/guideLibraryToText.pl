#!/usr/bin/perl

use Cwd;

my $inText = "GeCKO_lib_unique.txt";
my $outFile = "deduplicated_GeCKO_lib.txt";

open OUT, ">$outFile" or die $!;
open IN, "$inText" or die $!;

while (<IN>){
  my $inline = $_;
  $inline =~ s/\s+$//;
  my @inputArray = split /\t/, $inline;
  next if $inputArray[1] =~ 1;
  my $seq = $inputArray[2];
  my @infoArray = split /\|/, $inputArray[3];
  print "$infoArray[1]\t$seq\t$inputArray[2]\n";
  exit;
}
close OUT;
close IN;

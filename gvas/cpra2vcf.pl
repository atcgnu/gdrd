#! /usr/bin/perl

my ($cpra, $output) = @ARGV;

open OT, "> $output" or die $!;

print OT "##fileformat=VCFv4.0\n";
print OT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

my ($chr, $pos, $ref, $alt) = (split /[:>]/, $cpra);
#my ($ref, $alt) = (split /\>/, $ra);
$chr =~ s/chr//g;
print OT "chr$chr\t$pos\t.\t$ref\t$alt\t.\t.\t.\n";
close VCF;

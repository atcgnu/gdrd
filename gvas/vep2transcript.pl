#!/usr/bin/env perl

#***************************************************************************************************
# FileName: 
# Creator: Chen Y.L. <shenyulan@genomics.cn>
# Create Time: Thu Jun 26 17:50:59 CST 2014

# Description:
# CopyRight: 
# vision: 0.1
# ModifyList:
#
#  Revision: 0.1.9
#  Modifier: Chen Y. L. <shenyulan@genomics.cn>
#  ModifyTime: 20170911
#  ModifyReason:
#               fixed: a bug reported by wangjue@genomics.cn, variants in the same gene do not annotated with same canonical transcript (refer to https://atcgnu.myjetbrains.com/youtrack/issue/TR-3).
#***************************************************************************************************

use File::Basename;
use threads;
use DBI;
use FindBin qw($Bin);
use Getopt::Std;
use lib ("$FindBin::RealBin");


my $usage=<<usage;
        Usage: perl $0 <IN|vcf> <IN|config> <IN|sample.db> <IN| individual dir> <OT|output file>
        Example:
usage

die($usage) unless @ARGV >0;

my (%samples, @header, @sampleNames, %sites, %features);

my ($family_vcf, $out_result) = @ARGV;
my $result_dir = dirname $out_result;
my $result_file = basename $out_result;

my $md5vcf = `md5sum $family_vcf | /bin/awk '{print \$1}'`;
chomp $md5vcf;

print "Running the screening analyses...\n";
&screen('HEADER', $result_dir, $result_file);

for my $ssig_chr ('M',1..22,'X','Y')
#for my $ssig_chr (1..22,'X','Y')
#for my $ssig_chr (21)
{
    my $pid=fork();
    die "Cannot fork: $!" if (! defined $pid);
    if($pid ==0)
    {
        &screen($ssig_chr, $result_dir, "$ssig_chr.$md5vcf");
        exit;
    }
    else
    {
        push (@child_pids, $pid);
    }
}

for my $pid (@child_pids)
{
   waitpid($pid, 0);
}

for my $ssig_chr ('M', 1..22,'X', 'Y')
#for my $ssig_chr (1..22,'X', 'Y')
#for my $ssig_chr (21)
{
        `sort -t "      " -k2 -n $result_dir/$ssig_chr.$md5vcf.tsv >> ${out_result}.tsv;rm $result_dir/$ssig_chr.$md5vcf.tsv` if -e "$result_dir/$ssig_chr.$md5vcf.tsv";
}
print "The screening analyses are complete.\n";

sub screen {
my ($sig_chr, $out_dir, $out_file) = @_;

open VCF, "$family_vcf" or die $!;
open OT, ">$out_dir/$out_file.tsv" or die $!;

my (%ssite, %gcount, %gselected,  %sgene, %genes, %trs);

my $vep_key;
while(<VCF>){
        chomp;
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|AA_MAF|EA_MAF|EXON|INTRON|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|DISTANCE|STRAND|CLIN_SIG|CANONICAL|SYMBOL|SYMBOL_SOURCE|SIFT|PolyPhen|GMAF|BIOTYPE|ENSP|DOMAINS|CCDS|AFR_MAF|AMR_MAF|ASN_MAF|EUR_MAF|PUBMED">
        $vep_key = $1 if /##INFO=\<ID=CSQ,Number=\.,Type=String,Description="Consequence type as predicted by VEP\. Format: (\S+)"\>/;
        next if /^##/;
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  huangxinyuan7_control   guchengxiu1_case,##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|AA_MAF|EA_MAF|EXON|INTRON|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|DISTANCE|STRAND|CLIN_SIG|CANONICAL|SYMBOL|SYMBOL_SOURCE|SIFT|PolyPhen|GMAF|BIOTYPE|ENSP|DOMAINS|CCDS|AFR_MAF|AMR_MAF|ASN_MAF|EUR_MAF|PUBMED">
        print OT  "#CHR\tPOS\tREF\tAllele\tFeature\tBIOTYPE\tIs_Canonical\tHGNC(SYMBOL)\tConsequence\tIMPACT\tHGVSc\tHGVSp\tEXON\tINTRON\tSIFT\tPolyPhen\n" if /^#/ and $sig_chr eq 'HEADER';

    next if /^#/;


    my ($chr, $pos, $rs, $ref, $alt, $qua, $filter, $info, $format, @sample_infos) = (split /\t/);
        my $alt_r = $1 if $info =~ /ALT=(\S+?);/;
    $chr =~ s/chr//;
        next unless $chr eq $sig_chr;

    my @ks = (split /:/, $format);
        my $vep_info;
        map{$vep_info = $1 if /CSQ=(\S+)/;}(split /;/, $info);
#       $vep_info =~s/CSQ=//;
#       my $vep_key="Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|AA_MAF|EA_MAF|EXON|INTRON|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|DISTANCE|STRAND|CLIN_SIG|CANONICAL|SYMBOL|SYMBOL_SOURCE|SIFT|PolyPhen|GMAF|BIOTYPE|ENSP|DOMAINS|CCDS|HGVSc|HGVSp|AFR_MAF|AMR_MAF|ASN_MAF|EUR_MAF|PUBMED";
        my @vep_keys = (split /\|/,$vep_key);
#print "$vep_key\n";
        my (%feature_imps, %Canonical_transcripts, %symbols, %engenes);
        my ($s_consequence, $s_imp, $s_imp_no) = ('', '', 100);
        foreach my $vep_value (split /,/,$vep_info){
                my (@vep_values) = (split /\|/, $vep_value);
                foreach my $vep_k (0..$#vep_keys){
                        $sites{$vep_value}{$vep_keys[$vep_k]} = $vep_values[$vep_k] if $vep_values[$vep_k];
                }

                foreach my $consequence (split /\&/, $sites{$vep_value}{Consequence}){ # upstream_gene_variant
                        my ($imp, $imp_no);

                        ($imp, $imp_no) = (split /:/,$imps{$consequence});      # 'upstream_gene_variant' => "LOW:23",

                        if($s_imp_no > $imp_no){
                            $s_consequence = $consequence;
                            $s_imp = $imp;
                            $s_imp_no = $imp_no;
                        }

                        map{$sites{$vep_value}{"$_"} = '.' unless $sites{$vep_value}{"$_"} =~ /\S+/}('CLIN_SIG', $consequence, 'Consequence', 'CANONICAL', 'SYMBOL', 'Feature', 'BIOTYPE', 'HGVSc', 'HGVSp', 'EXON', 'INTRON', 'DOMAINS', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'SIFT', 'PolyPhen');

$sites{$vep_value}{HGVSc} =~ s/%3D/=/g;
$sites{$vep_value}{HGVSp} =~ s/%3D/=/g;
                        print OT "$sites{$vep_value}{SYMBOL}\t$sites{$vep_value}{Feature}\t$sites{$vep_value}{BIOTYPE}\t$sites{$vep_value}{CANONICAL}\t$sites{$vep_value}{Consequence}\t$imps{$s_consequence}\t$sites{$vep_value}{HGVSc}\t$sites{$vep_value}{HGVSp}\t'$sites{$vep_value}{EXON}'\t'$sites{$vep_value}{INTRON}'\t$sites{$vep_value}{SIFT}\t$sites{$vep_value}{PolyPhen}\n" if $sites{$vep_value}{SYMBOL} ne '.';
#                       print OT "$chr\t$pos\t$ref\t$alt\t$sites{$vep_value}{SYMBOL}\t$sites{$vep_value}{Feature}\t$sites{$vep_value}{BIOTYPE}\t$sites{$vep_value}{CANONICAL}\t$sites{$vep_value}{Consequence}\t$imps{$s_consequence}\t$sites{$vep_value}{HGVSc}\t$sites{$vep_value}{HGVSp}\t'$sites{$vep_value}{EXON}'\t'$sites{$vep_value}{INTRON}'\t$sites{$vep_value}{SIFT}\t$sites{$vep_value}{PolyPhen}\n" if $sites{$vep_value}{SYMBOL} ne '.';

                }
        }

    }
}
close VCF;

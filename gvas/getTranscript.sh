#! /bin/bash
cpra=$1
out=$2

export PATH=/data/gdrd/phase2/bin:$PATH
export PERL5LIB=/data/gdrd/phase2/software/perl524/lib:/data/gdrd/phase2/software/vcftools_0.1.12b/lib/perl5/site_perl

perl="/data/gdrd/phase2/software/perl524/bin/perl"
cpra2vcf="/data/gdrd/phase2/script/cpra2vcf.pl"
vep="/data/gdrd/phase2/software/ensembl-tools-release-77/scripts/variant_effect_predictor/variant_effect_predictor.pl"
vep_cache="/data/gdrd/phase2/data/cache"
screening="/data/gdrd/phase2/script/gvas.pl"

echo "transfer cpar to vcf"

$perl $cpra2vcf $cpra tmp.vcf
echo "VEP annotated"
$perl $vep -i tmp.vcf -o tmp.VEP.vcf --offline --dir_cache $vep_cache --vcf --force_overwrite --quiet --fork 6 --hgvs --assembly GRCh37 --everything
echo "screening..."
$perl $screening tmp.VEP.vcf $out

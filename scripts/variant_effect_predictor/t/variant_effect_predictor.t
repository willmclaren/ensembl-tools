use strict;
use warnings;
BEGIN { $| = 1;
	use Test::More;
	#plan tests => 6;
}

use FindBin qw($Bin);
use Data::Dumper;

#use Bio::EnsEMBL::Test::TestUtils;

# setup variables
my ($input, $output, $full_output, $expected);
my $tmpfile = "$$\_test_vep_input";

# configure script
my $script = $Bin.'/../variant_effect_predictor.pl';
my $perl   = '/usr/bin/env perl';
my $inc    = '-I ~/Variation/modules/';
my $ver    = 78;
my $ass    = 'GRCh38';
my $sp     = 'homo_sapiens';
my $cmd    = "$perl $inc $script -force -offline -dir_cache $Bin\/cache -i $tmpfile -o stdout -db $ver -assembly $ass -species $sp";

# unzip fasta
`gzip -d $Bin\/cache/$sp/$ver\_$ass/test.fa.gz` if(-e "$Bin\/cache/$sp/$ver\_$ass/test.fa.gz");


## BASIC TEST
#############

$output = `$cmd --help`;
ok($output =~ /ENSEMBL VARIANT EFFECT PREDICTOR/, "help message");


## HEADERS FOR OUTPUT FORMATS
#############################

# setup input
$input = '21 25587758 rs116645811 G A . . .';
input($input);

# vcf
$output = (grep {/^\#/} (split "\n", `$cmd -vcf`))[-1];
$expected = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
ok($output eq $expected, "output header - vcf") or diag("Expected\n$expected\n\nGot\n$output");

# default format
$full_output = `$cmd`;
$output = (grep {/^\#/} (split "\n", $full_output))[-1];
$expected =
  "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\t".
  "Feature_type\tConsequence\tcDNA_position\tCDS_position\t".
  "Protein_position\tAmino_acids\tCodons\tExisting_variation\tExtra";
  
ok($output eq $expected, "output header - default") or diag("Expected\n$expected\n\nGot\n$output");


## CONSEQUENCE TYPES
####################

# all consequence types should be covered by API tests
# just test a few to be sure

# missense
$output = (grep {/missense/} (split "\n", $full_output))[0];
$expected =
  "ENSG00000154719\tENST00000307301\tTranscript\tmissense_variant\t".
  "1043\t1001\t334\tT/M\taCg/aTg";

ok($output =~ /$expected/, "consequence type - missense") or diag("Expected\n$expected\n\nGot\n$output");

# upstream
$output = (grep {/upstream/} (split "\n", $full_output))[0];
$expected =
  "ENSG00000260583\tENST00000567517\tTranscript\tupstream_gene_variant\t".
  "-\t-\t-\t-\t-\t-\tDISTANCE=4432";

ok($output =~ /$expected/, "consequence type - upstream") or diag("Expected\n$expected\n\nGot\n$output");

# intron
$output = (grep {/intron/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000352957\tTranscript\tintron_variant";

ok($output =~ /$expected/, "consequence type - intron_variant") or diag("Expected\n$expected\n\nGot\n$output");


## VARIANT TYPES
################

# do all in one run to save time
$input = qq{21 25587758 25587758 C/T - rev
21 25587758 25587758 G/- + del
21 25587759 25587758 -/T + ins
21 25587758 25587760 GTA/TCC + multi
21 25587758 25587760 GTA/TC + unbalanced
21 25585656 25607517 DEL + svd 
21 25585656 25607517 DUP + svp 
};
input($input);
$full_output = `$cmd`;

# reverse strand
$output = (grep {/rev.+missense/} (split "\n", $full_output))[0];
$expected =
  "ENSG00000154719\tENST00000307301\tTranscript\tmissense_variant\t".
  "1043\t1001\t334\tT/M\taCg/aTg";
ok($output =~ /$expected/, "variant type - reverse strand") or diag("Expected\n$expected\n\nGot\n$output");

# deletion
$output = (grep {/del.+frameshift/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\tframeshift_variant\t1043\t1001\t334";
ok($output =~ /$expected/, "variant type - deletion") or diag("Expected\n$expected\n\nGot\n$output");

# insertion
$output = (grep {/ins.+frameshift/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\tframeshift_variant\t1042-1043\t1000-1001\t334";
ok($output =~ /$expected/, "variant type - insertion") or diag("Expected\n$expected\n\nGot\n$output");

# multi-bp sub
$output = (grep {/multi.+missense/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\tmissense_variant\t1041-1043\t999-1001\t333-334\tTT/TE\tacTACg/acGGAg";
ok($output =~ /$expected/, "variant type - multi-bp sub") or diag("Expected\n$expected\n\nGot\n$output");

# unbalanced sub
$output = (grep {/unbalanced.+frameshift/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\tframeshift_variant\t1041-1043\t999-1001\t333-334";
ok($output =~ /$expected/, "variant type - unbalanced sub") or diag("Expected\n$expected\n\nGot\n$output");

# sv - deletion
$output = (grep {/svd.+ENST00000307301/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\ttranscript_ablation";
ok($output =~ /$expected/, "variant type - structural variation - deletion") or diag("Expected\n$expected\n\nGot\n$output");

# sv - duplication
$output = (grep {/svp.+ENST00000307301/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\ttranscript_amplification";
ok($output =~ /$expected/, "variant type - structural variation - duplication") or diag("Expected\n$expected\n\nGot\n$output");


## OPTIONS
##########

## pathogenicity
$input = qq{21 25606454 25606454 G/C +};
input($input);
$full_output = `$cmd --sift b --polyphen b`;

# sift
$output = (grep {/ENST00000307301/} (split "\n", $full_output))[0];
$expected = "deleterious";
ok($full_output =~ /$expected/, "sift") or diag("Expected\n$expected\n\nGot\n$output");

# polyphen
$expected = "probably_damaging";
ok($output =~ /$expected/, "polyphen") or diag("Expected\n$expected\n\nGot\n$output");

# humdiv
$full_output = `$cmd --polyphen b --humdiv`;
$expected = 'probably_damaging\(0.998\)';
ok($full_output =~ /$expected/, "polyphen - humdiv") or diag("Expected\n$expected\n\nGot\n$full_output");


## hgvs
$full_output = `$cmd --hgvs`;
$output = (grep {/ENST00000419219/} (split "\n", $full_output))[0];
$expected = 'HGVSc=ENST00000419219.1:c.275C>G';
ok($output =~ /$expected/, "HGVSc") or diag("Expected\n$expected\n\nGot\n$output");
$expected = 'HGVSp=ENSP00000404426.1:p.Ala92Gly';
ok($output =~ /$expected/, "HGVSp") or diag("Expected\n$expected\n\nGot\n$output");


## regulation
$input = qq{21 25487468 25487468 A/T +};
input($input);
$full_output = `$cmd --regulatory`;

# reg feat
$output = (grep {/ENSR00000612061/} (split "\n", $full_output))[0];
$expected = 'regulatory_region_variant';
ok($output =~ /$expected/, "regfeat") or diag("Expected\n$expected\n\nGot\n$output");

# motif
$output = (grep {/MotifFeature/} (split "\n", $full_output))[0];
my %tmp_hash = split /\;|\=/, (split /\t/, $output)[-1];
$expected = {
  MOTIF_POS => 11,
  MOTIF_NAME => 'Name/Accession_association_EBF1:MA0154.2',
  HIGH_INF_POS => 'N',
  MOTIF_SCORE_CHANGE => -0.022,
  STRAND => 1
};
is_deeply(\%tmp_hash, $expected, "motif");

# cell type
$full_output = `$cmd --cell_type HUVEC`;
$output = (grep {/ENSR00000612061/} (split "\n", $full_output))[0];
$expected = 'promoter_flanking_region';
ok($output =~ /$expected/, "cell type") or diag("Expected\n$expected\n\nGot\n$output");


## colocated stuff
$input = qq{21 25584436 25584436 C/T +};
input($input);
$output = `$cmd --check_existing --gmaf --maf_1kg --pubmed`;

# colocated ID
$expected = 'rs2282471';
ok($output =~ /$expected/, "colocated ID") or diag("Expected\n$expected\n\nGot\n$output");

# gmaf
$expected = 'GMAF=T:0.1538';
ok($output =~ /$expected/, "GMAF") or diag("Expected\n$expected\n\nGot\n$output");

# population freqs
$expected = 'AFR_MAF=T:0.04';
ok($output =~ /$expected/, "AFR MAF") or diag("Expected\n$expected\n\nGot\n$output");
$expected = 'AMR_MAF=T:0.15';
ok($output =~ /$expected/, "AMR MAF") or diag("Expected\n$expected\n\nGot\n$output");
$expected = 'ASN_MAF=T:0.26';
ok($output =~ /$expected/, "ASN MAF") or diag("Expected\n$expected\n\nGot\n$output");
$expected = 'EUR_MAF=T:0.15';
ok($output =~ /$expected/, "EUR MAF") or diag("Expected\n$expected\n\nGot\n$output");

# pubmed
$expected = 'PUBMED=22272099';
ok($output =~ /$expected/, "pubmed") or diag("Expected\n$expected\n\nGot\n$output");

# somatic
$input = qq{21 25585742 25585742 G/A +};
input($input);
$output = `$cmd --check_existing`;
ok($output =~ /COSM162567.*SOMATIC\=1/, "somatic") or diag("Expected\n$expected\n\nGot\n$output");


## external IDs etc
$input = '21 25587758 rs116645811 G A . . .';
input($input);
$full_output = `$cmd --ccds --canonical --protein --uniprot --symbol --biotype --tsl --xref_refseq --numbers --domains`;
$output = (grep {/ENST00000352957/} (split "\n", $full_output))[0];

# gene symbol
$expected = 'SYMBOL=MRPL39';
ok($output =~ /$expected/, "gene symbol") or diag("Expected\n$expected\n\nGot\n$output");

# hgnc id
$expected = 'HGNC_ID=HGNC:14027';
ok($output =~ /$expected/, "HGNC ID") or diag("Expected\n$expected\n\nGot\n$output");

# biotype
$expected = 'BIOTYPE=protein_coding';
ok($output =~ /$expected/, "biotype") or diag("Expected\n$expected\n\nGot\n$output");

# CCDS
$expected = 'CCDS=CCDS13573.1';
ok($output =~ /$expected/, "CCDS") or diag("Expected\n$expected\n\nGot\n$output");

# protein
$expected = 'ENSP=ENSP00000284967';
ok($output =~ /$expected/, "protein ID") or diag("Expected\n$expected\n\nGot\n$output");

# swissprot
$expected = 'SWISSPROT=Q9NYK5';
ok($output =~ /$expected/, "SWISSPROT") or diag("Expected\n$expected\n\nGot\n$output");

# uniparc
$expected = 'UNIPARC=UPI00001AEE66';
ok($output =~ /$expected/, "UniParc") or diag("Expected\n$expected\n\nGot\n$output");

# refseq xref
$expected = 'RefSeq=NM_017446.3';
ok($output =~ /$expected/, "RefSeq xref") or diag("Expected\n$expected\n\nGot\n$output");

# TSL
$expected = 'TSL=1';
ok($output =~ /$expected/, "TSL") or diag("Expected\n$expected\n\nGot\n$output");

# numbers
$output = (grep {/ENST00000307301/} (split "\n", $full_output))[0];
$expected = 'EXON=10/11';
ok($output =~ /$expected/, "exon number") or diag("Expected\n$expected\n\nGot\n$output");

# protein domains
$expected = 'DOMAINS=Low_complexity_\(Seg\):Seg';
ok($output =~ /$expected/, "protein domains") or diag("Expected\n$expected\n\nGot\n$output");


## pick-type options
$full_output = `$cmd --pick`;
my @lines = grep {!/^\#/} (split "\n", $full_output);

ok(scalar @lines == 1, "pick - one line");
ok($lines[0] =~ /ENST00000307301/, "pick - correct transcript");

# per gene
$full_output = `$cmd --per_gene`;
@lines = grep {!/^\#/} (split "\n", $full_output);

ok(scalar @lines == 2, "per_gene");

# pick order
$input = '21 25716535 25716535 A/G +';
input($input);

my $o1 = `$cmd --pick --pick_order rank,length | grep -v '#'`;
my $o2 = `$cmd --pick | grep -v '#'`;

ok($o1 ne $o2 && $o1 =~ /ENST00000480456/ && $o2 =~ /ENST00000400532/, "pick order") or diag("--pick_order rank,length: $o1\ndefault order: $o2\n");

## INPUTS THAT CAUSE ERROR

# HGVS error
# 6 32191658 . TAGCAGCAGCAGC TAGCAGCAGCAGCAGC,T,TAGCAGCAGC 

#




## FINISHED
###########

unlink($tmpfile) if -e $tmpfile;
done_testing();


## SUBROUTINES
##############

# writes string $input to file in $tmpfile
sub input {
  my $data = shift;
  open OUT, ">$tmpfile" or die "ERROR: Could not write to temporary file $tmpfile\n";
  print OUT $data;
  close OUT;
}
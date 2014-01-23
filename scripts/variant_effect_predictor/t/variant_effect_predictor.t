use strict;
use warnings;
use Data::Dumper;
BEGIN { $| = 1;
	use Test::More;
	#plan tests => 6;
}

use FindBin qw($Bin);

#use Bio::EnsEMBL::Test::TestUtils;

# setup variables
my ($input, $output, $full_output, $expected);
my $tmpfile = "$$\_test_vep_input";

# configure script
my $script = $Bin.'/../variant_effect_predictor.pl';
my $perl   = '/usr/bin/env perl';
my $inc    = '-I ~/Variation/modules/';
my $cmd    = "$perl $inc $script -force -offline -dir_cache cache -i $tmpfile -o stdout -db 73";


## BASIC TEST
#############

$output = `$cmd --help`;
ok($output =~ /ENSEMBL VARIANT EFFECT PREDICTOR/, "help message");


## HEADERS FOR OUTPUT FORMATS
#############################

# setup input
$input = '21 26960070 rs116645811 G A . . .';
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

# gvf has no header, so no need to check


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
$input = qq{21 26960070 26960070 C/T - rev
21 26960070 26960070 G/- + del
21 26960071 26960070 -/T + ins
21 26960070 26960072 GTA/TCC + multi
21 26960070 26960072 GTA/TC + unbalanced
21 26957968 26979829 DEL + svd
21 26957968 26979829 DUP + svp
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
$expected = "ENSG00000154719\tENST00000307301\tTranscript\tframeshift_variant,feature_truncation\t1043\t1001\t334";
ok($output =~ /$expected/, "variant type - deletion") or diag("Expected\n$expected\n\nGot\n$output");

# insertion
$output = (grep {/ins.+frameshift/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\tframeshift_variant,feature_elongation\t1042-1043\t1000-1001\t334";
ok($output =~ /$expected/, "variant type - insertion") or diag("Expected\n$expected\n\nGot\n$output");

# multi-bp sub
$output = (grep {/multi.+missense/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\tmissense_variant\t1041-1043\t999-1001\t333-334\tTT/TE\tacTACg/acGGAg";
ok($output =~ /$expected/, "variant type - multi-bp sub") or diag("Expected\n$expected\n\nGot\n$output");

# unbalanced sub
$output = (grep {/unbalanced.+frameshift/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\tframeshift_variant,feature_truncation\t1041-1043\t999-1001\t333-334";
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
$input = qq{21 26978766 26978766 G/C +};
input($input);
$full_output = `$cmd --sift b --polyphen b`;

# sift
$output = (grep {/ENST00000307301/} (split "\n", $full_output))[0];
$expected = "deleterious";
ok($output =~ /$expected/, "sift") or diag("Expected\n$expected\n\nGot\n$output");

# polyphen
$expected = "probably_damaging";
ok($output =~ /$expected/, "polyphen") or diag("Expected\n$expected\n\nGot\n$output");


## regulation
$input = qq{21 26768619 26768619 T/G +};
input($input);
$full_output = `$cmd --regulatory`;

# reg feat
$output = (grep {/ENSR00000684049/} (split "\n", $full_output))[0];
$expected = 'regulatory_region_variant';
ok($output =~ /$expected/, "regfeat") or diag("Expected\n$expected\n\nGot\n$output");

# motif
$output = (grep {/MotifFeature/} (split "\n", $full_output))[0];
my %tmp_hash = split /\;|\=/, (split /\t/, $output)[-1];
$expected = {
  MOTIF_POS => 6,
  MOTIF_NAME => 'Jaspar_Matrix_NFKB:MA0105.1',
  HIGH_INF_POS => 'N',
  MOTIF_SCORE_CHANGE => -0.032
};
is_deeply(\%tmp_hash, $expected, "motif");

# cell type
$full_output = `$cmd --cell_type CD4`;
$output = (grep {/ENSR00000684049/} (split "\n", $full_output))[0];
$expected = 'Promoter_Associated';
ok($output =~ /$expected/, "cell type") or diag("Expected\n$expected\n\nGot\n$output");


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
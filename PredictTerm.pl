use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $genome_in = "NA";
my $rdt_in = "NA";
my $rit_in = "NA";
my $pred_in = "NA";
my $rdt_rf = "$cd_path/RFModels/rf_rdt_te.rds";
my $rit_rf = "$cd_path/RFModels/rf_rit_te.rds";
my $retrain = 0;
my $pe = 0;
my $gt = 0;
my $s_cutoff_rdt = 0.65;
my $s_cutoff_rit = 0.65;
my $tt_disc = 0.05;
my $window_shift = 100;
my $ref_bed = "NA";
my $rev_com = 0;

my $srcdir = $cd_path;

my $ind = 0;
foreach(@ARGV) {
	if (@ARGV[$ind] eq '-rdt_train') {
		$rdt_in = @ARGV[$ind + 1];
		$retrain = 1;
	}
	
	if (@ARGV[$ind] eq '-rit_train') {
		$rit_in = @ARGV[$ind + 1];
		$retrain = 1;
	}
	
	if (@ARGV[$ind] eq '-genome') {
		$genome_in = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-pred') {
		$pred_in = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rf_rdt') {
		$rdt_rf = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rf_rit') {
		$rit_rf = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-out') {
		$outdir = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-gt') {
		$gt = 1;
		$retrain = 0;
	}
	
	if (@ARGV[$ind] eq '-pe') {
		$pe = 1;
		$retrain = 0;
	}
	
	if (@ARGV[$ind] eq '-s_cutoff_rdt') {
		$s_cutoff_rdt = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-s_cutoff_rit') {
		$s_cutoff_rit = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-tt_disc') {
		$tt_disc  = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-window_shift') {
		$window_shift = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-ref_bed') {
		$ref_bed = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rc') {
		$rev_com = 1;
	}
	
	if (@ARGV[$ind] eq '-out') {
		$outdir = @ARGV[$ind + 1];
	}
	
	$ind++;
}

if (! -e $pred_in) {
	print "Input file $pred_in does not exit, exit.\n"; 
	exit;
}

if (-d $outdir) {
	system("rm -rf $outdir");
}
mkdir($outdir);

system("sh $cd_path/MakeWindows/ToWindows.sh $pred_in $outdir/MakeWindows $window_shift $rev_com > /dev/null 2> /dev/null");

if ($retrain > 0 and (! -e $rdt_in)) {
	$rdt_in = "$cd_path/test_set_ds/test_rdt_ds.fna";
}

if ($retrain > 0 and (! -e $rit_in)) {
	$rit_in = "$cd_path/test_set_ds/test_rit_ds.fna";
}

if ($pe > 0) {
	$rdt_rf = "$cd_path/RFModels/rf_rdt_te.rds";
	$rit_rf = "$cd_path/RFModels/rf_rit_te.rds";
}

if ($gt > 0 and $pe < 1) {
	$rdt_rf = "$cd_path/RFModels/rf_rdt_gt.rds";
	$rit_rf = "$cd_path/RFModels/rf_rit_gt.rds";
}

if (! -e $rdt_rf) {
	$rdt_rf = "$cd_path/RFModels/rf_rdt_te.rds";
}

if (! -e $rit_rf) {
	$rit_rf = "$cd_path/RFModels/rf_rit_te.rds";
}

if (-e $genome_in) {
	$genome_in = "$cd_path/test_set_ds/GCF_000750555.fna";
}

my $input = "$outdir/MakeWindows/input_seq.fna";
if ($retrain > 0) {
	system("perl $cd_path/TermClassify.pl -rdt_train $rdt_in -rit_train $rit_in -pred $input -out $outdir/TermClassifyMainProcess -genome $genome_in -mer 1 > /dev/null 2> /dev/null");
} else {
	system("perl $cd_path/TermClassify.pl -rf_rdt $rdt_rf -rf_rit $rit_rf -pred $input -out $outdir/TermClassifyMainProcess -genome $genome_in -mer 1 > /dev/null 2> /dev/null");
}

my $raw_csv = "$outdir/TermClassifyMainProcess/Pred/TermClassify/prob_rdt_rit_combined.csv";
my $bed_in = "$outdir/MakeWindows/input_seq.window.bed";
my $result_out = "$outdir/result";

if ((-e $raw_csv) and (-e $bed_in)) {
	system("Rscript $cd_path/PostProcess/TTCSVToBed_simple.r $raw_csv $bed_in $result_out $s_cutoff_rdt $s_cutoff_rit $tt_disc $ref_bed $rev_com > /dev/null 2> /dev/null");
} elsif (-e $raw_csv) {
	system("Rscript $cd_path/PostProcess/TTCSVToBed_simple.r $raw_csv $bed_in $result_out $s_cutoff_rdt $s_cutoff_rit $tt_disc $bed_in 0 > /dev/null 2> /dev/null");
}

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}
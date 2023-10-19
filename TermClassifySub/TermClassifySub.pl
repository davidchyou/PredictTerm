use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);
$cd_path_up = get_path($cd_path);

my $fa_in = "";
my $info_in = "NA";
my $rnie_gff = "NA";
my $outdir = "";
my $rf_rdt = "";
my $rf_rit = "";
my $te = 75;
my $bkgd = "NA";
my $feat_name_rdt = "rho_dependent";
my $feat_name_rit = "Independent";
my $user_data = "NA";

my $ind = 0;
foreach(@ARGV) {

	if (@ARGV[$ind] eq '-in') {
		$fa_in = @ARGV[$ind + 1];
		if (! (-e $fa_in)) {
			die "cannot open file: " . $fa_in . "\n";
		}
	}
	
	if (@ARGV[$ind] eq '-out') {
		$outdir = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-lookup') {
		$info_in = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rnie_gff') {
		$rnie_gff = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rf_rdt') {
		$rf_rdt = @ARGV[$ind + 1];
		if (! (-e $rf_rdt)) {
			die "cannot open file: " . $rf_rdt . "\n";
		}
	}
	
	if (@ARGV[$ind] eq '-rf_rit') {
		$rf_rit = @ARGV[$ind + 1];
		if (! (-e $rf_rit)) {
			die "cannot open file: " . $rf_rit . "\n";
		}
	}
	
	if (@ARGV[$ind] eq '-feat_name_rdt') {
		$feat_name_rdt = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-feat_name_rit') {
		$feat_name_rit = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-split') {
		$te = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-bkgd') {
		$bkgd = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-user_data') {
		$user_data = @ARGV[$ind + 1];
	}
	
	$ind++;
}

if (-d $outdir) {
	system("rm -rf $outdir");
}
mkdir($outdir);

if (-e $bkgd) {
	system("perl $cd_path_up/RUTSite/RunRUTSite.pl -in $fa_in -out $outdir/RDT -lookup $info_in -rnie_gff $rnie_gff -rf $rf_rdt -split $te -bkgd $bkgd -user_data $user_data");
	system("perl $cd_path_up/RUTSite/RunRUTSite.pl -in $fa_in -out $outdir/RIT -lookup $info_in -rnie_gff $rnie_gff -rf $rf_rit -split $te -bkgd $bkgd -user_data $user_data");
	system("Rscript $cd_path/TermClassifyMergeData.r $info_in $outdir/RDT/combined_prediction.csv $outdir/RIT/combined_prediction.csv $outdir/prob_rdt_rit_combined.csv $feat_name_rdt $feat_name_rit");
} else {
	system("perl $cd_path_up/RUTSite/RunRUTSite.pl -in $fa_in -out $outdir/RDT -lookup $info_in -rnie_gff $rnie_gff -rf $rf_rdt -split $te -user_data $user_data");
	system("perl $cd_path_up/RUTSite/RunRUTSite.pl -in $fa_in -out $outdir/RIT -lookup $info_in -rnie_gff $rnie_gff -rf $rf_rit -split $te -user_data $user_data");
	system("Rscript $cd_path/TermClassifyMergeData.r $info_in $outdir/RDT/combined_prediction.csv $outdir/RIT/combined_prediction.csv $outdir/prob_rdt_rit_combined.csv $feat_name_rdt $feat_name_rit");
}

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}
use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $fa_in = "";
my $lookup = "NA";
my $rnie_gff = "NA";
my $rf_in = "NA";
my $outdir = "";
my $train = 0;
my $data_only = 0;
my $te = 75;
my $bkgd = "NA";
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
		$lookup = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rnie_gff') {
		$rnie_gff = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rf') {
		$rf_in = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-train') {
		$train = 1;
	}
	
	if (@ARGV[$ind] eq '-data_only') {
		$data_only = 1;
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

system("Rscript $cd_path/RUTSite.r $fa_in $lookup $rnie_gff $rf_in $outdir $train $data_only $te $bkgd $user_data");

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}
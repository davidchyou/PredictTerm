use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $fa_in = "";
my $outdir = "";
my $split = 0;
my $downstream = 0;

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
	
	if (@ARGV[$ind] eq '-split') {
		$split = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-downstream') {
		$downstream = 1;
	}
	
	$ind++;
}

if (-d $outdir) {
	system("rm -rf $outdir");
}
mkdir($outdir);

my @toks = split("/", $fa_in);
my $ntoks = scalar(@toks);
my $fname = @toks[$ntoks - 1];

my $fa_work = $fa_in;
if ($split > 0) {
	system("Rscript $cd_path/Subseq.r $fa_work $outdir/$fname $split $downstream");
	$fa_work = "$outdir/$fname";
}

system("python $cd_path/DeSalvoRDT_predict.py $fa_work $outdir/$fname.orig.csv");
system("Rscript $cd_path/RTPPostProcess.r $outdir/$fname.orig.csv $fa_in $outdir/$fname.rtp.csv");

if (-e "$outdir/$fname") {
	unlink("$outdir/$fname");
}

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}
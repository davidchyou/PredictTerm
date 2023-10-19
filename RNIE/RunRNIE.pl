use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $fa_in = "";
my $outdir = "";
my $rnie_model_dir = "$cd_path/RNIE-master/models";
my $rnie_model = "gene.cm";
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
	
	if (@ARGV[$ind] eq '-rnie_md') {
		$rnie_model_dir = @ARGV[$ind + 1];
		if (! (-d $rnie_model_dir)) {
			die "cannot open RNIE CM-model directory: " . $rnie_model_dir . "\n";
		}
	}
	
	if (@ARGV[$ind] eq '-rnie_model') {
		$rnie_model =  @ARGV[$ind + 1];
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
my $fa_work = "$outdir/$fname.work.fa";
my $csv_out = "$outdir/$fname.rnie.csv";

my $curdir = getcwd();

if ($split > 0) {
	system("Rscript $cd_path/Subseq.r $fa_in $fa_work $split $downstream");
} else {
	system("cp $fa_in $fa_work");
}

chdir($outdir);
if (-e "$rnie_model_dir/$rnie_model") {
	system("cp $rnie_model_dir/$rnie_model gene.cm");
	#system("cmconvert -o gene.cm gene_orig.cm");
	system("perl $cd_path/RNIE-master/rnie.pl -f $fname.work.fa -gff -th 15 -o --toponly --gene -m gene.cm");
} else {
	system("cp $rnie_model gene.cm");
	#system("cmconvert -o gene.cm gene_orig.cm");
	system("perl $cd_path/RNIE-master/rnie.pl -f $fname.work.fa -gff -th 15 -o --toponly --gene -m gene.cm");
}

system("rm -rf $fname.work.fa");
chdir($curdir);

my @gffpaths = glob("$outdir/*rnie.gff");
my $gffpath = @gffpaths[0];
system("Rscript $cd_path/RNIEPostProcess.r $gffpath $fa_in $csv_out");

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}
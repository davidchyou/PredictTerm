use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $fa_in = "";
my $outdir = "";
my $genome_in = "";
my $train = 0;
my $rnd_seq = 0;
my $mer = 1;

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
	
	if (@ARGV[$ind] eq '-genome') {
		$genome_in = @ARGV[$ind + 1];
		#if (! (-e $genome_in)) {
		#	die "cannot open file: " . $genome_in . "\n";
		#}
	}
	
	if (@ARGV[$ind] eq '-train') {
		$train = 1;
	}
	
	if (@ARGV[$ind] eq '-rnd_seq') {
		$rnd_seq = 1;
	}
	
	if (@ARGV[$ind] eq '-mer') {
		$mer = @ARGV[$ind + 1];
	}
	
	$ind++;
}

if (-d $outdir) {
	system("rm -rf $outdir");
}
mkdir($outdir);

system("Rscript $cd_path/ProcessSeqs.r $fa_in $outdir $genome_in $train $rnd_seq $mer");

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}
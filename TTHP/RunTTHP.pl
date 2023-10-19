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

my $fa_work = "$outdir/$fname.work.fa";
my $csv_out = "$outdir/$fname.tt.csv";

my $curdir = getcwd();
system("cp $cd_path/transtermhp.py $outdir");
system("cp $cd_path/expterm.dat $outdir");

if ($split > 0) {
	system("Rscript $cd_path/Subseq.r $fa_in $fa_work $split $downstream");
} else {
	system("cp $fa_in $fa_work");
}

my $flag = chdir($outdir);
$flag = get_dummy_gff("$fname.work.fa", "dummy_gene.gff");

system("python ./transtermhp.py ./expterm.dat $fname.work.fa ./dummy_gene.gff tmp.gff 2>/dev/null");
system("grep terminator tmp.gff > $fname.tt.gff");

unlink("transtermhp.py");
unlink("expterm.dat");
unlink("tmp.coords");
unlink("tmp.fasta");
unlink("tmp.gff");
unlink("dummy_gene.gff");
unlink("$fname.work.fa");
$flag = chdir("$curdir");

my @gffpaths = glob("$outdir/*.tt.gff");
my $gffpath = @gffpaths[0];
system("Rscript $cd_path/TTHPPostProcess.r $gffpath $fa_in $csv_out");

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}

sub get_dummy_gff {
	my ($fa_in, $gff_out) = @_;
	
	open(IN, $fa_in);
	my @contents = <IN>;
	close(IN);
	
	open(GFF, ">>$gff_out");
	for (my $i = 0; $i < scalar(@contents); $i += 2) {
	
		my $head = @contents[$i];
		$head =~ s/\n//g;
		$head =~ s/^>//g;
	
		my $seq = @contents[$i + 1];
		$seq =~ s/\n//g;
		
		my $chrom_id = $head;
		my $len = length($seq);
		my $fgene2_start = $len - 1;
		my $fgene2_end = $len;
		
		print GFF "$chrom_id\tMYAPP\tgene\t1\t2\t.\t+\t1\tID=fakegene1\n";
		print GFF "$chrom_id\tMYAPP\tgene\t$fgene2_start\t$fgene2_end\t.\t+\t1\tID=fakegene2\n";
	}
	close(GFF);
	return 1;
}
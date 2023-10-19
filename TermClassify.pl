use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $genome_in = "NA";
my $genome_in_rdt = "NA";
my $genome_in_rit = "NA";
my $genome_in_pred = "NA";
my $gff_in = "NA";
my $gff_in_rdt = "NA";
my $gff_in_rit = "NA";
my $gff_in_pred = "NA";
my $rdt_in = "NA";
my $rit_in = "NA";
my $rdt_bed = "NA";
my $rit_bed = "NA";
my $pred_in = "NA";
my $pred_bed = "NA";
my $split_rdt = 75;
my $split_rit = 75;
my $split_pred = 75;
my $offset = 25;
my $outdir = "";
my $rdt_rf = "";
my $rit_rf = "";
my $rdt_rf_bkgd = "";
my $rit_rf_bkgd = "";
my $srcdir = $cd_path;
my $data_only = 0;
my $mer = 1;
my $bkgd_adj = 0;
my $run_rnie = 1;
my $non_seq_off_seq = 0;
my $run_other_tools = 0;
my $upto = 325;

my $ind = 0;
foreach(@ARGV) {
	if (@ARGV[$ind] eq '-rdt_train') {
		$rdt_in = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rit_train') {
		$rit_in = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rdt_pe_bed') {
		$rdt_bed = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rit_pe_bed') {
		$rit_bed = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-pred') {
		$pred_in = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-pred_pe_bed') {
		$pred_in = @ARGV[$ind + 1];
	}

	if (@ARGV[$ind] eq '-split_rdt') {
		$split_rdt = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-split_rit') {
		$split_rit = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-split_pred') {
		$split_pred = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rit_search_offset') {
		$offset = @ARGV[$ind + 1];
	}

	if (@ARGV[$ind] eq '-genome') {
		$genome_in = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-genome_rdt') {
		$genome_in_rdt = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-genome_rit') {
		$genome_in_rit = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-genome_pred') {
		$genome_in_pred = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-gff') {
		$gff_in = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-gff_rdt') {
		$gff_in_rdt = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-gff_rit') {
		$gff_in_rit = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-gff_pred') {
		$gff_in_pred = @ARGV[$ind + 1];
	}

	if (@ARGV[$ind] eq '-out') {
		$outdir = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-data_only') {
		$data_only = 1;
	}
	
	if (@ARGV[$ind] eq '-rf_rdt') {
		$rdt_rf = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rf_rit') {
		$rit_rf = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rf_rdt_bkgd') {
		$rdt_rf_bkgd = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-rf_rit_bkgd') {
		$rit_rf_bkgd = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-mer') {
		$mer = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-bkgd_adj') {
		$bkgd_adj = 1;
	}
	
	if (@ARGV[$ind] eq '-no_rnie') {
		$run_rnie = 0;
	}
	
	if (@ARGV[$ind] eq '-non_seq_off_seq') {
		$non_seq_off_seq = 1;
	}
	
	if (@ARGV[$ind] eq '-run_other_tools') {
		$run_other_tools = 1;
	}
	
	$ind++;
}

if (-d $outdir) {
	system("rm -rf $outdir");
}
mkdir($outdir);
 
if (-e $rdt_in) {
	my $rnie_split = $split_rdt + $offset;
	mkdir("$outdir/RDT_seqs");
	mkdir("$outdir/RDT");
	
	if (-e $genome_in_rdt) {
		$genome_in = $genome_in_rdt;
	}
	
	if (-e $gff_in_rdt) {
		$gff_in = $gff_in_rdt;
	}
	
	system("Rscript $srcdir/Util/HandlePEAndClass.r $rdt_in $rdt_bed $outdir/RDT_seqs/input_seq.fna $split_rdt rho-dependent 1");
	system("perl $srcdir/PreProcess/RunPreProcess.pl -in $outdir/RDT_seqs/input_seq.fna -out $outdir/RDT/PreProcess -genome $genome_in -train -mer $mer -rnd_seq");
	
	my $rnie_gff = "NA";
	if ($run_rnie > 0) {
		system("perl $srcdir/RNIE/RunRNIE.pl -in $outdir/RDT/PreProcess/input_seq.proc.fna -out $outdir/RDT/RNIE -split $rnie_split");
		$rnie_gff = "$outdir/RDT/RNIE/input_seq.proc.fna.work-geneMode-rnie.gff";
	}
	my $bkgd = "$outdir/RDT/PreProcess/input_seq.bkgd.fna";
	
	if ($run_other_tools > 0) {
		system("perl $srcdir/TTHP/RunTTHP.pl -in $outdir/RDT/PreProcess/input_seq.proc.fna -out $outdir/RDT/TTHP -split 0"); #-split $rnie_split
		system("perl $srcdir/RTP/RunRTP.pl -in $outdir/RDT/PreProcess/input_seq.proc.fna -out $outdir/RDT/RTP -split 0"); #-split $split_rdt -downstream
	}
	
	my $user_data = "NA";
	if ($non_seq_off_seq > 0 and -e $genome_in and -e $gff_in) {
		system("sh $srcdir/ROCO/ComputeRO.sh $outdir/RDT/PreProcess/input_seq.proc.fna $genome_in $gff_in $outdir/RDT/ROCO RDT $upto");
		$user_data = "$outdir/RDT/ROCO/ROCO.csv";
	}
	
	if ($bkgd_adj > 0) {
		if ($data_only < 1) {
			system("perl $srcdir/RUTSite/RunRUTSite.pl -in $outdir/RDT/PreProcess/input_seq.proc.fna -out $outdir/RDT/RFModel -lookup $outdir/RDT/PreProcess/input_seq.lookup.csv -rnie_gff $rnie_gff -train -split $split_rdt -bkgd $bkgd -user_data $user_data");
			$rdt_rf_bkgd = "$outdir/RDT/RFModel/rf_model.rds";
		} else {
			system("perl $srcdir/RUTSite/RunRUTSite.pl -in $outdir/RDT/PreProcess/input_seq.proc.fna -out $outdir/RDT/RFData -lookup $outdir/RDT/PreProcess/input_seq.lookup.csv -rnie_gff $rnie_gff -train -split $split_rdt -bkgd $bkgd -data_only -user_data $user_data");
		}
	} else {
		if ($data_only < 1) {
			system("perl $srcdir/RUTSite/RunRUTSite.pl -in $outdir/RDT/PreProcess/input_seq.proc.fna -out $outdir/RDT/RFModel -lookup $outdir/RDT/PreProcess/input_seq.lookup.csv -rnie_gff $rnie_gff -train -split $split_rdt -user_data $user_data");
			$rdt_rf = "$outdir/RDT/RFModel/rf_model.rds";
		} else {
			system("perl $srcdir/RUTSite/RunRUTSite.pl -in $outdir/RDT/PreProcess/input_seq.proc.fna -out $outdir/RDT/RFData -lookup $outdir/RDT/PreProcess/input_seq.lookup.csv -rnie_gff $rnie_gff -train -split $split_rdt -data_only -user_data $user_data");
		}
	}
	
	if ($run_other_tools > 0) {
		system("Rscript $srcdir/Util/AddData.r $outdir/RDT/RFModel/combined_validation.csv $outdir/RDT/TTHP/input_seq.proc.fna.tt.csv $outdir/RDT/RTP/input_seq.proc.fna.rtp.csv $outdir/RDT/RFModel/combined_validation_appended.csv"); 
	}
}

if (-e $rit_in) {
	my $rnie_split = $split_rit + $offset;
	mkdir("$outdir/RIT_seqs");
	mkdir("$outdir/RIT");
	
	if (-e $genome_in_rit) {
		$genome_in = $genome_in_rit;
	}
	
	if (-e $gff_in_rit) {
		$gff_in = $gff_in_rit;
	}
	
	system("Rscript $srcdir/Util/HandlePEAndClass.r $rit_in $rit_bed $outdir/RIT_seqs/input_seq.fna $split_rit Independent 1");
	system("perl $srcdir/PreProcess/RunPreProcess.pl -in $outdir/RIT_seqs/input_seq.fna -out $outdir/RIT/PreProcess -genome $genome_in -train -mer $mer -rnd_seq");
	
	my $rnie_gff = "NA";
	if ($run_rnie > 0) {
		system("perl $srcdir/RNIE/RunRNIE.pl -in $outdir/RIT/PreProcess/input_seq.proc.fna -out $outdir/RIT/RNIE -split $rnie_split");
		$rnie_gff = "$outdir/RIT/RNIE/input_seq.proc.fna.work-geneMode-rnie.gff";
	}
	my $bkgd = "$outdir/RIT/PreProcess/input_seq.bkgd.fna";
	
	if ($run_other_tools > 0) {
		system("perl $srcdir/TTHP/RunTTHP.pl -in $outdir/RIT/PreProcess/input_seq.proc.fna -out $outdir/RIT/TTHP -split 0"); #-split $rnie_split
		system("perl $srcdir/RTP/RunRTP.pl -in $outdir/RIT/PreProcess/input_seq.proc.fna -out $outdir/RIT/RTP -split 0"); #-split $split_rdt -downstream
	}
	
	my $user_data = "NA";
	if ($non_seq_off_seq > 0 and -e $genome_in and -e $gff_in) {
		system("sh $srcdir/ROCO/ComputeRO.sh $outdir/RIT/PreProcess/input_seq.proc.fna $genome_in $gff_in $outdir/RIT/ROCO IT $upto");
		$user_data = "$outdir/RIT/ROCO/ROCO.csv";
	}
	
	if ($bkgd_adj > 0) {
		if ($data_only < 1) {
			system("perl $srcdir/RUTSite/RunRUTSite.pl -in $outdir/RIT/PreProcess/input_seq.proc.fna -out $outdir/RIT/RFModel -lookup $outdir/RIT/PreProcess/input_seq.lookup.csv -rnie_gff $rnie_gff -train -split $split_rdt -bkgd $bkgd -user_data $user_data");
			$rit_rf_bkgd = "$outdir/RIT/RFModel/rf_model.rds";
		} else {
			system("perl $srcdir/RUTSite/RunRUTSite.pl -in $outdir/RIT/PreProcess/input_seq.proc.fna -out $outdir/RIT/RFData -lookup $outdir/RIT/PreProcess/input_seq.lookup.csv -rnie_gff $rnie_gff -train -split $split_rdt -bkgd $bkgd -data_only -user_data $user_data");
		}
	} else {
		if ($data_only < 1) {
			system("perl $srcdir/RUTSite/RunRUTSite.pl -in $outdir/RIT/PreProcess/input_seq.proc.fna -out $outdir/RIT/RFModel -lookup $outdir/RIT/PreProcess/input_seq.lookup.csv -rnie_gff $rnie_gff -train -split $split_rdt -user_data $user_data");
			$rit_rf = "$outdir/RIT/RFModel/rf_model.rds";
		} else {
			system("perl $srcdir/RUTSite/RunRUTSite.pl -in $outdir/RIT/PreProcess/input_seq.proc.fna -out $outdir/RIT/RFData -lookup $outdir/RIT/PreProcess/input_seq.lookup.csv -rnie_gff $rnie_gff -train -split $split_rdt -data_only -user_data $user_data");
		}
	}
	
	if ($run_other_tools > 0) {
		system("Rscript $srcdir/Util/AddData.r $outdir/RIT/RFModel/combined_validation.csv $outdir/RIT/TTHP/input_seq.proc.fna.tt.csv $outdir/RIT/RTP/input_seq.proc.fna.rtp.csv $outdir/RIT/RFModel/combined_validation_appended.csv"); 
	}
}

if (-e $pred_in) {
	my $rnie_split = $split_pred + $offset;
	mkdir("$outdir/Pred_seqs");
	mkdir("$outdir/Pred");
	
	if (-e $genome_in_pred) {
		$genome_in = $genome_in_pred;
	}
	
	if (-e $gff_in_pred) {
		$gff_in = $gff_in_pred;
	}
	
	system("Rscript $srcdir/Util/HandlePEAndClass.r $pred_in $pred_bed $outdir/Pred_seqs/input_seq.fna $split_pred NA 0");
	system("perl $srcdir/PreProcess/RunPreProcess.pl -in $outdir/Pred_seqs/input_seq.fna -out $outdir/Pred/PreProcess -genome $genome_in -mer $mer -rnd_seq");
	
	my $rnie_gff = "NA";
	if ($run_rnie > 0) {
		system("perl $srcdir/RNIE/RunRNIE.pl -in $outdir/Pred/PreProcess/input_seq.proc.fna -out $outdir/Pred/RNIE -split $rnie_split");
		$rnie_gff = "$outdir/Pred/RNIE/input_seq.proc.fna.work-geneMode-rnie.gff";
	}
	my $bkgd = "$outdir/Pred/PreProcess/input_seq.bkgd.fna";
	
	if ($run_other_tools > 0) {
		system("perl $srcdir/TTHP/RunTTHP.pl -in $outdir/Pred/PreProcess/input_seq.proc.fna -out $outdir/Pred/TTHP -split 0"); # -split $rnie_split
		system("perl $srcdir/RTP/RunRTP.pl -in $outdir/Pred/PreProcess/input_seq.proc.fna -out $outdir/Pred/RTP -split 0"); #$split_rdt -downstream
	}
	
	my $user_data = "NA";
	if ($non_seq_off_seq > 0 and -e $genome_in and -e $gff_in) {
		system("sh $srcdir/ROCO/ComputeRO.sh $outdir/Pred/PreProcess/input_seq.proc.fna $genome_in $gff_in $outdir/Pred/ROCO GENOME $upto");
		$user_data = "$outdir/Pred/ROCO/ROCO.csv";
	}
	
	if ($data_only > 0) {
		if ($bkgd_adj > 0) {
			system("perl $srcdir/RUTSite/RunRUTSite.pl -in $outdir/Pred/PreProcess/input_seq.proc.fna -out $outdir/Pred/PredData -lookup $outdir/Pred/PreProcess/input_seq.lookup.csv -rnie_gff $rnie_gff -split $split_rdt -bkgd $bkgd -data_only -user_data $user_data");
		} else {
			system("perl $srcdir/RUTSite/RunRUTSite.pl -in $outdir/Pred/PreProcess/input_seq.proc.fna -out $outdir/Pred/PredData -lookup $outdir/Pred/PreProcess/input_seq.lookup.csv -rnie_gff $rnie_gff -split $split_rdt -data_only -user_data $user_data");
		}
		exit;
	}
	
	if ($bkgd_adj > 0) { 
		$rdt_rf = $rdt_rf_bkgd;
		$rit_rf = $rit_rf_bkgd;
	}
	
	if ((-e $rdt_rf) && (-e $rit_rf)) {
		if ($bkgd_adj > 0) {
			system("perl $srcdir/TermClassifySub/TermClassifySub.pl -in $outdir/Pred/PreProcess/input_seq.proc.fna -out $outdir/Pred/TermClassify -lookup $outdir/Pred/PreProcess/input_seq.lookup.csv -rnie_gff $rnie_gff -rf_rdt $rdt_rf -rf_rit $rit_rf -split $split_pred -bkgd $bkgd -user_data $user_data");
		} else {
			system("perl $srcdir/TermClassifySub/TermClassifySub.pl -in $outdir/Pred/PreProcess/input_seq.proc.fna -out $outdir/Pred/TermClassify -lookup $outdir/Pred/PreProcess/input_seq.lookup.csv -rnie_gff $rnie_gff -rf_rdt $rdt_rf -rf_rit $rit_rf -split $split_pred -user_data $user_data");
		}
	}
}

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}
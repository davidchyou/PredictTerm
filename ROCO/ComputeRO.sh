SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
CURDIR=`pwd -P`
INFILE=$1 #"/mnt/SSD/brownlab/RDT_RIT_2nd_wave_paper/RDT_RIT_Output_paper_2/test_TermClassify_1mer/RDT/PreProcess/input_seq.proc.fna"
GENOME=$2 #"/mnt/SSD/brownlab/RDT_RIT_2nd_wave_paper/test_data_paper/GCF_000750555.fna"
GFF=$3 #"/mnt/SSD/brownlab/RDT_RIT_2nd_wave_paper/test_data_paper/GCF_000750555.cds.gff"
OUTDIR=$4 #"ROCO"
MODE=$5
SEQLEN=$6 #325

rm -rf $OUTDIR
mkdir $OUTDIR

$SCRIPTPATH/blastn -query $INFILE -subject $GENOME -outfmt '6 std sstrand' -max_target_seqs 1 2>/dev/null | \
awk '$4=='$SEQLEN | sort -k1,1 -u > $OUTDIR/input_seq_blast.bed
awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\n", $2, $9, $9+1,$1, $12, $13=="plus"?"+":"-"}' $OUTDIR/input_seq_blast.bed | \
sort -k1,1 -k2,2n > $OUTDIR/input_seq_blast_start.bed
awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\n", $2, $10, $10+1,$1, $12, $13=="plus"?"+":"-"}' $OUTDIR/input_seq_blast.bed | \
sort -k1,1 -k2,2n > $OUTDIR/input_seq_blast_end.bed
awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\n", $2, ($10>$9)?$9:$10, ($10>$9)?$10:$9, $1, $12, $13=="plus"?"+":"-"}' $OUTDIR/input_seq_blast.bed | \
sort -k1,1 -k2,2n > $OUTDIR/input_seq_blast_intv.bed
awk -F'\t' '$3=="CDS"' $GFF | \
awk -F'\t' '{printf "%s\t%s\t%s\t%s_%s\t%s\t%s\n", $1, $4, $4+1, $1, NR, 1, $7}' | sort -k1,1 -k2,2n > $OUTDIR/all_cds_start.bed
awk -F'\t' '$3=="CDS"' $GFF | \
awk -F'\t' '{printf "%s\t%s\t%s\t%s_%s\t%s\t%s\n", $1, $5, $5+1, $1, NR, 1, $7}' | sort -k1,1 -k2,2n> $OUTDIR/all_cds_end.bed
awk -F'\t' '$3=="CDS"' $GFF | \
awk -F'\t' '{printf "%s\t%s\t%s\t%s_%s\t%s\t%s\n", $1, $4, $5, $1, NR, 1, $7}' | sort -k1,1 -k2,2n> $OUTDIR/all_cds.bed
$SCRIPTPATH/bedtools closest -a $OUTDIR/input_seq_blast_start.bed -b $OUTDIR/all_cds_end.bed | sort -k4,4 -u | \
awk '$6=="-" {printf "%s\t%s\t%s\t%s\n", $4, $6, $10, $12}' > $OUTDIR/rev_term.txt
$SCRIPTPATH/bedtools closest -a $OUTDIR/input_seq_blast_end.bed -b $OUTDIR/all_cds_start.bed | sort -k4,4 -u | \
awk '$6=="+" {printf "%s\t%s\t%s\t%s\n", $4, $6, $10, $12}' > $OUTDIR/fwd_term.txt
cat $OUTDIR/fwd_term.txt $OUTDIR/rev_term.txt | sort -k1,1 -u | awk '{printf "%s\t%s\t%s\n", $1, $3=="+"?1:0, $4=="+"?1:0}' | \
awk '{printf "%s\t%s\n", $1, $2-$3}' | awk '{printf "%s\t%s\n", $1, $2==0?"0":"1"}' > $OUTDIR/conv_div.txt
Rscript $SCRIPTPATH/CDSOverlap.r $INFILE $GENOME $GFF $OUTDIR/input_seq_blast.bed $OUTDIR/conv_div.txt $MODE $OUTDIR/ROCO.csv

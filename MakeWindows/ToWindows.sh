SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
CURDIR=`pwd -P`
INFILE=$1 #"../test_set_ds/test_all_ds.fna"
OUTDIR=$2 #"out"
SHIFT=$3
REVCOM=$4
if [[ $REVCOM -lt 1 ]]
then
	rm -rf $OUTDIR; mkdir $OUTDIR; cp $INFILE $OUTDIR/input_seq.orig.fna
	perl -pe 's/>(.*)/ >\1\t/g; s/\n//g; s/\s>/\n>/g' $OUTDIR/input_seq.orig.fna | grep -v '^$' | perl -F'\t' -lane '{$s1=$F[0];$s2="NNNNNNNNNN".$F[1];print "$s1\n$s2"}' > $OUTDIR/input_seq.orig1.fna
	perl -pe 's/>(.*)/ >\1\t/g; s/\n//g; s/\s>/\n>/g' $OUTDIR/input_seq.orig1.fna | grep -v '^$' | perl -F'\t' -lane '{$s1=$F[0];$s2=$F[1];$s1=~s/^>//g;$ll=length($s2);print "$s1\t10\t$ll\t$s1\t1\t+"}' > $OUTDIR/input_seq.bed
	$SCRIPTPATH/bedtools makewindows -b $OUTDIR/input_seq.bed -w 325 -s 100 -i srcwinnum | awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, 1, "+"}' > $OUTDIR/input_seq.window1.bed
	$SCRIPTPATH/bedtools getfasta -fi $OUTDIR/input_seq.orig1.fna -bed $OUTDIR/input_seq.window1.bed -nameOnly -s | perl -pe 's/\([+-]\)$//g' > $OUTDIR/input_seq.fna
	awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1, ($2-10<1)?1:($2-10), $3-10, $4, 1, "+"}' $OUTDIR/input_seq.window1.bed > $OUTDIR/input_seq.window.bed
	rm -rf $OUTDIR/input_seq.orig1.fna; rm -rf $OUTDIR/input_seq.window1.bed rm -rf $OUTDIR/input_seq.orig1.fna.fai
	#| perl -pe 's/(\S+)::\S+/\1/g'
else
	rm -rf $OUTDIR; mkdir $OUTDIR; cp $INFILE $OUTDIR/input_seq.orig.fna
	perl -pe 's/>(.*)/ >\1\t/g; s/\n//g; s/\s>/\n>/g' $OUTDIR/input_seq.orig.fna | grep -v '^$' | perl -F'\t' -lane '{$s1=$F[0];$s2="NNNNNNNNNN".$F[1];print "$s1\n$s2"}' > $OUTDIR/input_seq.orig1.fna
	perl -pe 's/>(.*)/ >\1\t/g; s/\n//g; s/\s>/\n>/g' $OUTDIR/input_seq.orig1.fna | grep -v '^$' | perl -F'\t' -lane '{$s1=$F[0];$s2=$F[1];$s1=~s/^>//g;$ll=length($s2);print "$s1\t10\t$ll\t$s1\t1\t-"}' > $OUTDIR/input_seq.bed
	$SCRIPTPATH/bedtools makewindows -b $OUTDIR/input_seq.bed -w 325 -s 100 -i srcwinnum | awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, 1, "-"}' > $OUTDIR/input_seq.window1.bed
	$SCRIPTPATH/bedtools getfasta -fi $OUTDIR/input_seq.orig1.fna -bed $OUTDIR/input_seq.window1.bed -nameOnly -s | perl -pe 's/\([+-]\)$//g' > $OUTDIR/input_seq.fna
	awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1, ($2-10<1)?1:($2-10), $3-10, $4, 1, "-"}' $OUTDIR/input_seq.window1.bed > $OUTDIR/input_seq.window.bed
	rm -rf $OUTDIR/input_seq.orig1.fna; rm -rf $OUTDIR/input_seq.window1.bed rm -rf $OUTDIR/input_seq.orig1.fna.fai
	#| perl -pe 's/(\S+)::\S+/\1/g'
fi

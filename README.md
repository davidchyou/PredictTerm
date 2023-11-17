# PredictTerm

**What it does?**

PredictTerm predicts Rho-dependent (RDT) and Intrinsic Terminators (IT) from a bacterial genome sequence. It is based on the set of RDT features published by Nadiras et al., combining with IT prediction scores from the IT-prediction software RNIE, together with features around mRNA processed-ends that distinguishes IT and RDT from non-terminators, and between the two. Two random forests are trained to score for RDT and IT based on the aforementioned parameters, the RDT and IT scores were then used to classify terminators into IT, RDT, both IT and RDT (IT+RDT) and unclassified. Unclassified sequences can be non-terminators, or novel terminators when sequences are transcription terminator extracts. This is a repository for a publication in preparation.

Highlights: Predicts both main types of transcription terminators in E coli and scores results.

**Dependencies**

All executables are provided, including the version of Infernal used by the IT-predictor RNIE as well as other RNIE dependencies. R libraries need to be installed separately using the "install.packages" function. The "Biostrings" package is part of Bioconductor which involves the installation of "BiocManager" (install.packages("BiocManager")) if not previously installed.

- RNIE (provided)
- Bedtools (provided)
- R library stringr (install.packages("stringr"))
- R library zoo (install.packages("zoo"))
- R library dplyr (install.packages("dplyr"))
- R library Biostrings (BiocManager::install("Biostrings"))
- R library randomForest (install.packages("randomForest"))

**Basis uses**

When input is a set of genomic extracts in multi-FASTA format, specifying the input file path and the path of output directory will suffice.

        perl PredictTerm.pl -pred test_set_ju/NC_000913.tt.pe.fna -out NC_000913_out_pe

By default, PredictTerm uses a set of random forests that are trained on RDT and IT extracts around the processed ends (PE), spanning PE-74 to PE+250. However, PredictTerm will slide along the sequences in windows of 325nt, so exact position of PE is not essential. If sequences were gene-termini extracts, there is another set of random forests for these and PredictTerm will use them with the "-gt" argument.

        perl PredictTerm.pl -gt -pred test_set_ju/NC_000913.tt.gt.fna -out NC_000913_out_gt

By default, PredictTerm will report the position of terminator, the type, RDT and IT scores in BED6 format, relative to the input sequences. If input sequences are extracted from a genome and genomic positions of input sequences are available in BED6 format, PredictTerm will convert terminator positions into true genomic coordinates by providing it using the "-ref_bed" argument. 

        perl PredictTerm.pl -pred test_set_ju/NC_000913.tt.pe.fna -out NC_000913_out_pe -ref_bed test_set_ju/NC_000913.tt.pe.bed

**Scanning the whole genome**

Because PredictTerm can slide along sequences to predict terminators and call types, and report the results in BED6 format, it has the capacity to scan for terminators in a genome and call types. PredictTerm will digitize sequences longer than 325nt into windows of 325nt, and the argument "-window_shift" can be used to specify the step-size.

        perl PredictTerm.pl -window_shift 40 -pred test_set_ju/NC_000913.fna -out NC_000913_out_plus
        perl PredictTerm.pl -window_shift 40 -rc -pred test_set_ju/NC_000913.fna -out NC_000913_out_minus

Don't forget the reverse strand! (-rc)

Main input (detail below): fna file, [gff file of CDS positions]

Main output (detail below): Terminator predictions site, type, scores 

Advanced features you may want to change: Can retain on different genome, e.g. one with a different GC%.

What to do with the output: Predicts strong and weak terminators weak terminators may have regulatory roles. Understanding the role of Rho and RDT in termination.

What is similar? PredictTerm is unique in predict both IT and RDT is ease of use and accuracy.  Previous IT predictors: RNIE, TransTermHP, RDT predictor: RhoTermPredict.

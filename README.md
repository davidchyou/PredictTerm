# PredictTerm

**What it does**

PredictTerm predicts Rho-dependent (RDT) and Intrinsic Terminators (IT) from a bacterial genome sequence. It is based on some known RDT features (Nadiras et al.2019), and with IT prediction scores from the IT-prediction software RNIE (Gardner et al 2011), together with features around mRNA processed-ends that distinguishes IT and RDT from non-terminators, and between the two. Two random forests are trained to score for RDT and IT based on the aforementioned parameters, the RDT and IT scores were then used to classify terminators into IT, RDT, both IT and RDT (IT+RDT) and unclassified. Unclassified sequences can be non-terminators, or novel terminators when sequences are transcription terminator extracts. This is a repository for a publication in preparation.

Highlights: Predicts both main types of transcription terminators in E coli and scores results.

**Dependencies**

All executables are provided, including the older version of Infernal used by the IT-predictor RNIE as well as other RNIE dependencies. R libraries need to be installed separately using the "install.packages" function. The "Biostrings" package is part of Bioconductor which involves the installation of "BiocManager" (install.packages("BiocManager")) if not previously installed.

- RNIE (provided)
- Bedtools (provided)
- R library stringr (install.packages("stringr"))
- R library zoo (install.packages("zoo"))
- R library dplyr (install.packages("dplyr"))
- R library Biostrings (BiocManager::install("Biostrings"))
- R library randomForest (install.packages("randomForest"))

**Basic use**

When input is a set of genomic extracts in multi-FASTA format, specifying the input file path and the path of output directory will suffice.

        perl PredictTerm.pl -pred test_set_ju/NC_000913.tt.pe.fna -out NC_000913_out_pe

By default, PredictTerm uses a set of random forests that are trained on RDT and IT extracts around the processed ends (PE), spanning PE-74 to PE+250. However, PredictTerm will slide along the sequences in windows of 325nt, so exact position of PE is not essential. If sequences were gene-termini extracts, there is another set of random forests for these and PredictTerm will use them with the "-gt" argument.

        perl PredictTerm.pl -gt -pred test_set_ju/NC_000913.tt.gt.fna -out NC_000913_out_gt

By default, PredictTerm will report the position of terminator, the type, RDT and IT scores in BED6 format, relative to the input sequences. If input sequences are extracted from a genome and genomic positions of input sequences are available in BED6 format, PredictTerm will convert terminator positions into true genomic coordinates by providing it using the "-ref_bed" argument. 

        perl PredictTerm.pl -pred test_set_ju/NC_000913.tt.pe.fna -out NC_000913_out_pe -ref_bed test_set_ju/NC_000913.tt.pe.bed

**Scanning a whole genome**

Because PredictTerm can slide along sequences to predict terminators and call types, and report the results in BED6 format, it has the capacity to scan for terminators in a genome and call types. PredictTerm will digitize sequences longer than 325nt into windows of 325nt, and the argument "-window_shift" can be used to specify the step-size.

        perl PredictTerm.pl -window_shift 40 -rc -pred test_set_ju/NC_000913.fna -out NC_000913_out_plus
        perl PredictTerm.pl -window_shift 40 -rc -pred test_set_ju/NC_000913.fna -out NC_000913_out_minus

Don't forget the reverse strand! (-rc)

**Retraining and prediction**

Random Forests in PredictTerm were trained based on E coli terminators published by Dar and Sorek 2016. The application as a whole were tested on terminators in a different E coli strain reported by Ju et al. For non-E coli terminators, users can use the training layer to train a different set of RDT and IT random forests based on the training set given, and then predict terminators and call types using the trained random forests. The arugments "-rdt_train" and "-rit_train" specifies the RDT and IT training sets respectively. The argument "-genome" specifies the genome where the training sets were derived, PredictTerm will generated shuffled genomic extracts from that to create the negative set. If the "-genome" argument is not used, negative training set will be generated randomly where nucleotide compositions are all 0.25. Here is an example.

         perl PredictTerm.pl -rit_train test_set_ds/test_rit_ds.fna -rdt_train test_set_ds/test_rdt_ds.fna -genome test_set_ds/GCF_000750555.fna test_set_ds -pred test_set_ju/NC_000913.tt.pe.fna -out NC_000913_out_retrain -ref_bed test_set_ju/NC_000913.tt.pe.bed

Leave-half-out cross-validation will be performed during training. ROC-curves and AUC statistics will be generated for reporting. Retraining will not replace the default random forests (which were trained based on RDT and IT reported by Dar and Sorek).  

**What to do with the output?**

PredictTerm can be used to predict strong terminators and type them as IT, RDT, IT+RDT or unclassified. While "unclassified terminators" may be non-terminators, these can be novel terminators if within terminator regions or near gene termini. Terminators found within transcribed regions or operons may have regulatory roles. As an example, we can use PredictTerm to understand the role of Rho factor and RDT in termination. 

**Related Software**

PredictTerm is unique in the way that it predicts both IT and RDT and classifies them accurately. It is also easy of use. 
Here is a list of related software that we comapared in this deveoplment.

- RNIE: This IT predictor called by PredictTerm as a dependency to assist with conventional IT prediction. It uses Infernal and pretrained coveriation model of intrinsic terminators.
- TransTermHP: This is a pattern-based intrinsic terminator predictor. The accuracy is good but requires non-standard data files for CDS annotations.
- RhoTermPredict: This is a pattern-based RDT predictor, it scans for RDT in genomes or sequences efficiently with moderately-level of accuracy. The RDT pattern used by RhoTermPredict was also considered by Nadiras's RDT predictor, and random forests of PredictTerm uses these as part of the parameter set.
- Nadiras RDT classifier: The RDT classifier implemented by Nadiras et al. collects features and patterns relevant to RDT. Partial least-square discriminant analysis was then used to classify sequences into strong RDT, weak RDT and non-RDT with excellent reported accuracy. The Python script generating the RDT parameter is publicly available, but the classifier is private and requires paid software (SIMCA) to run.
- iTerm-PseKNC: This terminator predictor uses physicochemical and pseudonucleotide compositions to predict terminators, but does not call types. It is highly accurate on intrinsic terminators.
- BATTER (unpublished but available as an unreviewed preprint on BioRxiv): It is a "terminator or not" terminator predictor like iTerm-PseKNC. It is based on a language model of terminator stem-loops, including both RDT and IT. The preprint reported accuracy is excellent for in RDT and IT, and also terminators in genomes with biased nucleotide distributions. However, although BATTER predicts both RDT and IT as terminators, it does not distinguish between RDT and IT, nor handling the case where an IT stem-loop is upstream to a RDT site (IT+RDT). Terminator type calling is the competitive advantage of PredictTerm.

Users may consider using iTerm-PseKNC or BATTER to scan for terminators in genomes or sequences, and then apply PredictTerm on putative terminator regions to determine terminator types.

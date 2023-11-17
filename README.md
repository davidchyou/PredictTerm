# PredictTerm

**What it does?**
PredictTerm predicts Rho-dependent (RDT) and Intrinsic Terminators (IT) from a bacterial genome sequence. It is based on the set of RDT features published by Nadiras et al., combining with IT prediction scores from the IT-prediction software RNIE, together with features around mRNA processed-ends that distinguishes IT and RDT from non-terminators, and between the two. Two random forests are trained to score for RDT and IT based on the aforementioned parameters, the RDT and IT scores were then used to classify terminators into IT, RDT, both IT and RDT (IT+RDT) and unclassified. Unclassified sequences can be non-terminators, or novel terminators when sequences are transcription terminator extracts. This is a repository for a publication in preparation.

Highlights: Predicts both main types of transcription terminators in E coli and scores results.

**Dependencies**
All executables are provided, particularly the version of Infernal used by the IT-predictor RNIE.

How to get it and dependencies (details below): github

Example of use with test files: perl PredictTerm.pl -pred input_genome.fna -out input_genome.out 
perl PredictTerm.pl -pred NC_000913.fna -out NC_000913.out

Main input (detail below): fna file, [gff file of CDS positions]

Main output (detail below): Terminator predictions site, type, scores

Advanced features you may want to change: Can retain on different genome, e.g. one with a different GC%.

What to do with the output: Predicts strong and weak terminators weak terminators may have regulatory roles. Understanding the role of Rho and RDT in termination.

What is similar? PredictTerm is unique in predict both IT and RDT is ease of use and accuracy.  Previous IT predictors: RNIE, TransTermHP, RDT predictor: RhoTermPredict.

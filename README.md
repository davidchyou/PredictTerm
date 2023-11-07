# PredictTerm

Description: Predicts Rho-dependent (RDT) and Intrinsic Terminators (IT) from a bacterial genome sequence. This is a repository for a publication in preparation.

Highlights: Predicts both main types of transcription terminators in E coli, scores results.

How to get it and dependencies (details below): github

Example of use with test files: perl PredictTerm.pl -pred input_genome.fna -out input_genome.out 

Main input (detail below): fna file, [gff file of CDS positions]

Main output (detail below): Terminator predictions site, type, scores

Advanced features you may want to change: Can retain on different genome, e.g. one with a different GC%.

What to do with the output: Predicts strong and weak terminators weak terminators may have regulatory roles. Understanding the role of Rho and RDT in termination.

What is similar? PredictTerm is unique in predict both IT and RDT is ease of use and accuracy.  Previous IT predictors: RNIE, TransTermHP, RDT predictor: RhoTermPredict.

/* cmalign.c
 * SRE, Thu Jul 25 11:28:03 2002 [St. Louis]
 * SVN $Id: cmalign.c 3013 2009-10-28 14:09:25Z nawrockie $
 * 
 * Align sequences to a CM.
 * 
 *****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************
 */
#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"		/* general seq analysis library   */
#include "esl_alphabet.h"
#include "esl_getopts.h"		
#include "esl_mpi.h"
#include "esl_msa.h"
#include "esl_random.h"		
#include "esl_sq.h"		
#include "esl_sqio.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

#define ALGOPTS  "--cyk,--optacc,--viterbi,--sample" /* Exclusive choice for algorithm */
#define OUTALPHOPTS "--rna,--dna"                    /* Exclusive choice for output alphabet */
#define ACCOPTS  "--nonbanded,--hbanded,--qdb"       /* Exclusive choice for acceleration strategies */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-o",        eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "output the alignment to file <f>, not stdout", 1 },
  { "-l",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "align locally w.r.t. the model",         1 },
  { "-p",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,   "--small", "append posterior probabilities to alignment", 1 },
  { "-q",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "quiet; suppress banner and scores, print only the alignment", 1 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "output alnment in non-interleaved, 1 line/seq Stockholm format",  1 },
  { "--informat",eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL,        NULL, "specify the input file is in format <x>, not FASTA", 1 },
  { "--devhelp", eslARG_NONE,   NULL,  NULL, NULL,      NULL,      NULL,        NULL, "show list of undocumented developer options", 1 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    "--merge", "run as an MPI parallel program",                    1 },  
#endif
  /* --merge: merge two alignments and quit, don't do any new aligning */
  { "--merge",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "merge two alignments and exit", 8 },
  /* Algorithm options */
  { "--optacc",  eslARG_NONE,"default", NULL, NULL,     ALGOPTS,    NULL,"--small,--qdb", "align with the Holmes/Durbin optimal accuracy algorithm", 2 },
  { "--cyk",     eslARG_NONE,   FALSE,  NULL, NULL,     ALGOPTS,    NULL,        NULL, "align with the CYK algorithm", 2 },
  { "--sample",  eslARG_NONE,   FALSE,  NULL, NULL,     ALGOPTS,    NULL,"--small,--qdb", "sample alignment of each seq from posterior distribution", 2 },
  { "-s",        eslARG_INT,     NULL, NULL, "n>0",      NULL,"--sample",        NULL, "w/--sample, set random number generator seed to <n>",  2 },
  { "--viterbi", eslARG_NONE,   FALSE,  NULL, NULL,     ALGOPTS,    NULL,        "-p", "align to a CM Plan 9 HMM with the Viterbi algorithm",2 },
  { "--sub",     eslARG_NONE,   FALSE,  NULL, NULL,     NULL,       NULL,        "-l", "build sub CM for columns b/t HMM predicted start/end points", 2 },
  { "--small",   eslARG_NONE,   FALSE,  NULL, NULL,     NULL,"--cyk",     "--hbanded", "use divide and conquer (d&c) alignment algorithm", 2 },
  /* Banded alignment */
  { "--hbanded", eslARG_NONE, "default",  NULL, NULL,   NULL,     NULL,      ACCOPTS, "accelerate using CM plan 9 HMM derived bands", 3 },
  { "--nonbanded",eslARG_NONE,  FALSE, NULL, NULL,"--hbanded",    NULL,  "--hbanded", "do not use bands to accelerate aln algorithm", 3 },
  { "--tau",     eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,"--hbanded",       NULL, "set tail loss prob for --hbanded to <x>", 3 },
  { "--mxsize",  eslARG_REAL, "2048.0",NULL, "x>0.",     NULL,      NULL,"--small,--qdb", "set maximum allowable DP matrix size to <x> Mb", 3},
  /* Options that modify how the output alignment is created */
  { "--rna",     eslARG_NONE,"default",NULL, NULL,  OUTALPHOPTS,   NULL,        NULL, "output alignment as RNA sequence data", 4},
  { "--dna",     eslARG_NONE,   FALSE, NULL, NULL,  OUTALPHOPTS,   NULL,        NULL, "output alignment as DNA (not RNA) sequence data", 4},
  { "--matchonly",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        "-p", "include only match columns in output alignment", 4 },
  { "--resonly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "include only match columns with >= 1 residues in output aln", 4 },
  { "--fins",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "flush inserts left/right in output alignment", 4 },
  { "--onepost", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      "-p",        NULL, "with -p, only append single '0-9,*' char as posterior probability", 4 },
  /* Including a preset alignment */
  { "--withali", eslARG_INFILE, NULL,  NULL, NULL,      NULL,      NULL,  "--viterbi","incl. alignment in <f> (must be aln <cm file> was built from)", 5 },
  { "--withpknots",eslARG_NONE, NULL,  NULL, NULL,      NULL,"--withali",       NULL, "incl. structure (w/pknots) from <f> from --withali <f>", 5 },
  { "--rf",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--withali",       NULL, "--rf was originally used with cmbuild", 5 },
  { "--gapthresh",eslARG_REAL,  "0.5", NULL, "0<=x<=1", NULL,"--withali",       NULL, "--gapthresh <x> was originally used with cmbuild", 5 },
  /* Verbose output files */
  { "--tfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "dump individual sequence parsetrees to file <f>", 7 },

  /* All options below are developer options, only shown if --devhelp invoked */
  /* Developer options related to alignment algorithm */
  { "--inside",   eslARG_NONE,  FALSE, NULL, NULL,      ALGOPTS,   NULL,     ALGOPTS, "don't align; return scores from the Inside algorithm", 101 },
  { "--checkpost",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      "-p",        NULL, "check that posteriors are correctly calc'ed", 101 },
  { "--no-null3",eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "turn OFF the NULL3 post hoc additional null model", 101 },
  /* developer options related to banded alignment */
  { "--checkfb", eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       "-l", "check that HMM posteriors for bands were correctly calc'ed", 102},
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "use posterior sums during HMM band calculation (widens bands)", 102 },
  { "--qdb",     eslARG_NONE,   FALSE, NULL, NULL,"--hbanded",     NULL,     ACCOPTS, "use query dependent banded CYK alignment algorithm", 102 },
  { "--beta",    eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,   "--qdb",        NULL, "set tail loss prob for --qdb to <x>", 102 },
  { "--hsafe",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded","--viterbi,-p,--optacc", "realign (w/o bands) seqs with HMM banded CYK score < 0 bits", 102 },
  /* developer options related to output files and debugging */
  { "--regress", eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save regression test data to file <f>", 103 },
  { "--banddump",eslARG_INT,    "0",   NULL, "0<=n<=3", NULL,      NULL,        NULL, "set verbosity of band info print statements to <n>", 103 },
  { "--dlev",    eslARG_INT,    "0",   NULL, "0<=n<=3", NULL,      NULL,        NULL, "set verbosity of debugging print statements to <n>", 103 },
  { "--stall",   eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,        NULL, "arrest after start: for debugging MPI under gdb", 103 },  
  /* Developer options related to experiment local begin/end modes */
  { "--pebegin", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      "-l",  "--pbegin", "set all local begins as equiprobable", 104 },
  { "--pfend",   eslARG_REAL,   NULL,  NULL, "0<x<1",   NULL,      "-l",    "--pend", "set all local end probs to <x>", 104 },
  { "--pbegin",  eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      "-l",        NULL, "set aggregate local begin prob to <x>", 104 },
  { "--pend",    eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      "-l",        NULL, "set aggregate local end prob to <x>", 104 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  char         *cmfile;	        /* name of input CM file  */ 
  char         *sqfile;	        /* name of sequence file  */ 
  ESL_SQFILE   *sqfp;           /* open sequence input file stream */
  int           fmt;		/* format code for seqfile */
  ESL_ALPHABET *abc;		/* digital alphabet for the CM */
  int           ncm;            /* number cm we're on */

  int           do_mpi;		/* TRUE if we're doing MPI parallelization */
  int           nproc;		/* how many MPI processes, total */
  int           my_rank;	/* who am I, in 0..nproc-1 */
  int           do_stall;	/* TRUE to stall the program until gdb attaches */
  ESL_RANDOMNESS *r;            /* source of randomness, only created if --sample enabled */

  /* Masters only (i/o streams) */
  CMFILE       *cmfp;		/* open input CM file stream       */
  FILE         *ofp;		/* output file (default is stdout) */
  FILE         *tracefp;	/* optional output for parsetrees  */
  FILE         *regressfp;	/* optional output for regression test  */
  ESL_MSAFILE  *withalifp;	/* optional input alignment to include */
  ESL_MSA      *withmsa;	/* MSA from withalifp to include */
  char         *withss_cons;	/* ss_cons string from withmsa (before knot stripping) */
  Parsetree_t  *withali_mtr;	/* guide tree for MSA from withalifp */
  ESL_ALPHABET *withali_abc;	/* digital alphabet for reading withali MSA */
  ESL_ALPHABET *abc_out;	/* digital alphabet for output */
};

static char usage1[] = "[-options] <cmfile> <sequence file>";
static char usage2[] = "[-options] --merge <cmfile> <msafile1> <msafile2>";
static char banner[] = "align sequences to an RNA CM";

static int  init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
/* static int  init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf); */

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  serial_merge_alignments_only(const ESL_GETOPTS *go);
#ifdef HAVE_MPI
static int   mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int   mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif

static int process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln);
static int output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln);

static int initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int check_withali(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, ESL_MSA **ret_msa, Parsetree_t **ret_mtr);
static int include_withali(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr, int *ret_nseq, char *errbuf);
static int compare_cm_guide_trees(CM_t *cm1, CM_t *cm2);
static int make_aligned_string(char *aseq, char *gapstring, int alen, char *ss, char **ret_s);
static int add_withali_pknots(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_MSA *newmsa);

static int print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
static void print_cm_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command);
static int get_markup_from_msa(ESL_MSA *msa, char *errbuf, char ***ret_comment, int *ret_ncomment, char ***ret_gf_tag, char ***ret_gf, int *ret_ngf, char ***ret_gs_tag, char ****ret_gs, int *ret_ngs, char ***ret_gr_tag, char ****ret_gr, int *ret_ngr);
static int append_parsetrees(Parsetree_t **tr_to_append, int ntr_to_append, Parsetree_t ***ret_tr, int *ret_ntr, char *errbuf);
static int append_sequences(ESL_SQ **sq_to_append, int nsq_to_append, ESL_SQ ***ret_sq, int *ret_nsq, char *errbuf);
static int add_msa_markup(ESL_MSA *merged_msa, char *errbuf, int nseq, int seq_offset,
			  char **comment, int ncomment,
			  char **gf_tag, char **gf,  int ngf,
			  char **gs_tag, char ***gs, int ngs,
			  char **gr_tag, char ***gr, int ngr);
extern int gapize_string_to_fit_alignment(char *s, const char *aseq, int alen, char gapchar_to_add, char *aln_gapchars, char *errbuf, char **ret_news);
			    

/*
  static void print_stage_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg);
  static int print_align_options(const struct cfg_s *cfg, CM_t *cm);
*/

#ifdef HAVE_MPI
static int determine_nseq_per_worker(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, int *ret_nseq_worker);
static int add_worker_seqs_to_master(seqs_to_aln_t *master_seqs, seqs_to_aln_t *worker_seqs, int offset);
#endif

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("Memory error, stopwatch not created.\n");
  esl_stopwatch_Start(w);
  struct cfg_s     cfg;
  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

  /*********************************************** 
   * Parse command line
   ***********************************************/

  /* Process command line options.
   */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || 
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage1);
      esl_usage(stdout, argv[0], usage2);
      printf("\n  The --merge option merges the two alignments in <msafile1> and <msafile2>\n  created by previous runs of cmalign with <cmfile> into a single alignment.");
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "--devhelp") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage1);
      esl_usage(stdout, argv[0], usage2);
      puts("\n  The --merge option merges the two alignments in <msafile1> and <msafile2>\n  created by previous runs of cmalign with <cmfile> into a single alignment.");
      puts("\nwhere general options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\nalignment algorithm related options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nbanded dynamic programming acceleration options:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\noutput options:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\nmerge alignments in <msafile1> and <msafile2>:");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
      puts("\noptions for including a fixed alignment within output alignment:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\nverbose output files and debugging:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      puts("\nundocumented developer algorithm options:");
      esl_opt_DisplayHelp(stdout, go, 101, 2, 80);
      puts("\nundocumented developer banded alignment options:");
      esl_opt_DisplayHelp(stdout, go, 102, 2, 80);
      puts("\nundocumented developer verbose output/debugging options:");
      esl_opt_DisplayHelp(stdout, go, 103, 2, 80);
      puts("\nundocumented developer options for experimental local begin/end modes:");
      esl_opt_DisplayHelp(stdout, go, 104, 2, 80);
      exit(0);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage1);
      esl_usage(stdout, argv[0], usage2);
      puts("\n  The --merge option merges the two alignments in <msafile1> and <msafile2>\n  created by previous runs of cmalign with <cmfile> into a single alignment.");
      puts("\nwhere general options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\nalignment algorithm related options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nbanded dynamic programming acceleration options:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\noutput options:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\nmerge alignments in <msafile1> and <msafile2>:");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
      puts("\noptions for including a fixed alignment within output alignment:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\nverbose output files and debugging:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      exit(0);
    }
  if(((! esl_opt_GetBoolean(go, "--merge")) && (esl_opt_ArgNumber(go) != 2)) ||
     ((  esl_opt_GetBoolean(go, "--merge")) && (esl_opt_ArgNumber(go) != 3))) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage1);
      esl_usage(stdout, argv[0], usage2);
      puts("\n  The --merge option merges the two alignments in <msafile1> and <msafile2>\n  created by previous runs of cmalign with <cmfile> into a single alignment.");
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      printf("\nTo see more help on other available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  /* Check for incompatible option combinations I don't know how to disallow with esl_getopts */
  /* --small requires EITHER --nonbanded or --qdb */
  if ((esl_opt_GetBoolean(go, "--small")) && (! ((esl_opt_GetBoolean(go, "--nonbanded")) || (esl_opt_GetBoolean(go, "--qdb"))))) { 
    printf("Error parsing options, --small is only allowed in combination with --nonbanded or --qdb.\n");
    exit(1);
  }

  /* if --merge, merge the two alignments and exit, never create cfg structure we won't need it */
  if(esl_opt_GetBoolean(go, "--merge")) { 
    serial_merge_alignments_only(go);
    esl_stopwatch_Destroy(w);
    esl_getopts_Destroy(go);
    return 0;	
  }

  /* Initialize what we can in the config structure (without knowing the input alphabet yet).
   */
  cfg.cmfile     = esl_opt_GetArg(go, 1); 
  cfg.sqfile     = esl_opt_GetArg(go, 2); 
  cfg.sqfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  if   (esl_opt_IsDefault(go, "--informat")) cfg.fmt = eslSQFILE_UNKNOWN; /* autodetect sequence file format by default. */ 
  else { 
    cfg.fmt = esl_sqio_FormatCode(esl_opt_GetString(go, "--informat"));
    if(cfg.fmt == eslSQFILE_UNKNOWN) cm_Fail("Can't recognize sequence file format: %s. valid options are: fasta, embl, genbank, ddbj, uniprot, stockholm, or pfam\n", esl_opt_GetString(go, "--informat"));
  }
  cfg.abc        = NULL;	           /* created in init_master_cfg() in masters, or in mpi_worker() in workers */
  if      (esl_opt_GetBoolean(go, "--rna")) cfg.abc_out = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna")) cfg.abc_out = esl_alphabet_Create(eslDNA);
  else    cm_Fail("Can't determine output alphabet");
  cfg.cmfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.ofp        = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.tracefp    = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.regressfp  = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.withalifp  = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.withmsa    = NULL;	           /* filled in init_master_cfg() in masters, stays NULL for workers */
  cfg.withss_cons= NULL;	           /* filled in check_withali() in masters, stays NULL for workers */
  cfg.withali_mtr= NULL;	           /* filled in init_master_cfg() in masters, stays NULL for workers */
  cfg.withali_abc= NULL;	           /* created in init_master_cfg() in masters, stays NULL for workers */
  cfg.ncm        = 0;
  cfg.r          = NULL;	           /* created in init_master_cfg() for masters, mpi_worker() for workers*/

  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */
  cfg.do_stall   = esl_opt_GetBoolean(go, "--stall");


  /* This is our stall point, if we need to wait until we get a
   * debugger attached to this process for debugging (especially
   * useful for MPI):
   */
  while (cfg.do_stall); 

  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HAVE_MPI
  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      int              status;               /* easel status */
      char             errbuf[cmERRBUFSIZE]; /* for error messages in mpi_master() */
      cfg.do_mpi     = TRUE;

      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if(cfg.nproc == 1) cm_Fail("ERROR, MPI mode, but only 1 processor running...");

      if (cfg.my_rank > 0) { status = mpi_worker(go, &cfg); }
      else { 
	if(! esl_opt_GetBoolean(go, "-q")) cm_banner(stdout, argv[0], banner);
	status = mpi_master(go, &cfg, errbuf);
      }
      /* check status, if eslOK, we continue, else we exit. either way we call MPI_Finalize() */
      if(status == eslOK) { 
	esl_stopwatch_Stop(w);
	esl_stopwatch_MPIReduce(w, 0, MPI_COMM_WORLD);
	MPI_Finalize();
      }
      else { /* status != eslOK, master has error message in errbuf, worker does not */
	MPI_Finalize();
	if(cfg.my_rank == 0) cm_Fail(errbuf); /* master */
	else                 return 0;        /* worker */
      }
    }
  else
#endif /*HAVE_MPI*/
    {
      if(! esl_opt_GetBoolean(go, "-q")) cm_banner(stdout, argv[0], banner);
      serial_master(go, &cfg);
      esl_stopwatch_Stop(w);
    }
  /* Clean up the shared cfg. 
   */
  if (cfg.my_rank == 0) {
    if (! esl_opt_IsDefault(go, "-o")) { 
      printf("# Alignment saved in file %s.\n", esl_opt_GetString(go, "-o"));
      fclose(cfg.ofp); 
    }
    if (cfg.tracefp   != NULL) { 
      printf("# Parsetrees saved in file %s.\n", esl_opt_GetString(go, "--tfile"));
      fclose(cfg.tracefp);
    }
    if (cfg.regressfp   != NULL) {
      printf("# Regression data (alignment) saved in file %s.\n", esl_opt_GetString(go, "--regress"));
      fclose(cfg.regressfp);
    }
    if (cfg.cmfp      != NULL) CMFileClose(cfg.cmfp);
    if (cfg.sqfp      != NULL) esl_sqfile_Close(cfg.sqfp);
    if (cfg.withalifp != NULL) esl_msafile_Close(cfg.withalifp);
    if (cfg.withmsa   != NULL) esl_msa_Destroy(cfg.withmsa);
    if (cfg.withali_mtr != NULL) FreeParsetree(cfg.withali_mtr);
    if (cfg.withss_cons != NULL) free(cfg.withss_cons);
  }
  if (cfg.r         != NULL) esl_randomness_Destroy(cfg.r);
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.abc_out   != NULL) esl_alphabet_Destroy(cfg.abc_out);
  if (cfg.withali_abc != NULL) esl_alphabet_Destroy(cfg.withali_abc);
  if (cfg.my_rank == 0 && (! esl_opt_GetBoolean(go, "-q"))) { 
    printf("#\n");
    esl_stopwatch_Display(stdout, w, "# CPU time: ");
  }
  esl_getopts_Destroy(go);
  esl_stopwatch_Destroy(w);
  return 0;
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Already set:
 *    cfg->cmfile      - command line arg 1
 *    cfg->sqfile      - command line arg 2
 *    cfg->fmt         - format of output file
 * Sets: 
 *    cfg->sqfp        - open sequence file                
 *    cfg->ofp         - output file (stdout by default)
 *    cfg->cmfp        - open CM file                
 *    cfg->abc         - digital input alphabet
 *    cfg->tracefp     - optional output file
 *    cfg->regressfp   - optional output file
 *    cfg->withalifp   - optional input alignment file to include
 *    cfg->withmsa     - MSA from --withali file 
 *    cfg->withali_mtr - guide tree for MSA from --withali file 
 *    cfg->withali_abc - digital input alphabet for --withali file
 *    cfg->r           - source of randomness
 *                   
 * Errors in the MPI master here are considered to be "recoverable",
 * in the sense that we'll try to delay output of the error message
 * until we've cleanly shut down the worker processes. Therefore
 * errors return (code, errbuf) by the ESL_FAIL mech.
 */
static int
init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;
  int type;

  /* open input sequence file */
  status = esl_sqfile_Open(cfg->sqfile, cfg->fmt, NULL, &(cfg->sqfp));
  if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "File %s doesn't exist or is not readable\n", cfg->sqfile);
  else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of sequence file %s\n", cfg->sqfile);
  else if (status == eslEINVAL)  ESL_FAIL(status, errbuf, "Can’t autodetect stdin or .gz."); 
  else if (status != eslOK)      ESL_FAIL(status, errbuf, "Sequence file open failed with error %d\n", status);
  if(cfg->sqfp->format == eslMSAFILE_STOCKHOLM) ESL_FAIL(eslEFORMAT, errbuf, "cmalign doesn't support Stockholm alignment format. Please reformat to FASTA.\n");
  cfg->fmt = cfg->sqfp->format;

  /* Set the sqfile alphabet as RNA, if it's DNA we're fine. 
   * If it's not RNA nor DNA, we can't deal with it anyway,
   * so we're hardcoded to RNA.
   */
  cfg->abc = esl_alphabet_Create(eslRNA);
  if(cfg->abc == NULL) ESL_FAIL(status, errbuf, "Failed to create alphabet for sequence file");
  esl_sqfile_SetDigital(cfg->sqfp, cfg->abc);

  /* open CM file */
  if((cfg->cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
   ESL_FAIL(eslFAIL, errbuf, "Failed to open covariance model save file %s\n", cfg->cmfile);

  /* open output file */
  if (esl_opt_GetString(go, "-o") != NULL) {
    if ((cfg->ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
    } else cfg->ofp = stdout;

  /* seed master's RNG, this will only be used if --sample enabled, but we always initialize it for convenience (seeds always get sent to workers) */
  if (! esl_opt_IsDefault(go, "-s")) 
    cfg->r = esl_randomness_Create((long) esl_opt_GetInteger(go, "-s"));
  else cfg->r = esl_randomness_CreateTimeseeded();

  /* optionally, open trace file */
  if (esl_opt_GetString(go, "--tfile") != NULL) {
    if ((cfg->tracefp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
    }

  /* optionally, open regression file */
  if (esl_opt_GetString(go, "--regress") != NULL) {
    if ((cfg->regressfp = fopen(esl_opt_GetString(go, "--regress"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --regress output file %s\n", esl_opt_GetString(go, "--regress"));
    }

  /* optionally, open withali file for reading */
  if(esl_opt_GetString(go, "--withali") != NULL)
    {
      status = esl_msafile_Open(esl_opt_GetString(go, "--withali"), eslMSAFILE_UNKNOWN, NULL, &(cfg->withalifp));
      if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "--withali alignment file %s doesn't exist or is not readable\n", 
					      esl_opt_GetString(go, "--withali"));
      else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of --withali alignment %s\n", 
					      esl_opt_GetString(go, "--withali"));
      else if (status != eslOK)      ESL_FAIL(status, errbuf, "Alignment file open failed with error %d\n", status);
      /* Guess the withali alphabet, if it's ambiguous, guess RNA,
       * we'll treat RNA and DNA both as RNA internally.
       * We can't handle any other alphabets, so this is hardcoded. */
      status = esl_msafile_GuessAlphabet(cfg->withalifp, &type);
      if (status == eslEAMBIGUOUS)    type = eslRNA; /* guess it's RNA, we'll fail downstream with an error message if it's not */
      else if (status == eslEFORMAT)  ESL_FAIL(status, errbuf, "Alignment file parse failed: %s\n", cfg->withalifp->errbuf);
      else if (status == eslENODATA)  ESL_FAIL(status, errbuf, "Alignment file %s is empty\n", esl_opt_GetString(go, "--withali"));
      else if (status != eslOK)       ESL_FAIL(status, errbuf, "Failed to read alignment file %s\n", esl_opt_GetString(go, "--withali"));
      /* we can read DNA/RNA but internally we treat it as RNA */
      if(! (type == eslRNA || type == eslDNA))
	ESL_FAIL(eslEFORMAT, errbuf, "Alphabet is not DNA/RNA in %s\n", esl_opt_GetString(go, "--withali"));
      cfg->withali_abc = esl_alphabet_Create(eslRNA);
      if(cfg->withali_abc == NULL) ESL_FAIL(status, errbuf, "Failed to create alphabet for --withali");
      esl_msafile_SetDigital(cfg->withalifp, cfg->withali_abc);
    }

  if(cfg->r == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create master RNG.");

  return eslOK;
}

/* serial_master()
 * The serial version of cmalign.
 * Align each sequence to the CM.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;
  char     errbuf[cmERRBUFSIZE];
  CM_t     *cm;
  seqs_to_aln_t  *seqs_to_aln;  /* sequences to align, holds seqs, parsetrees, CP9 traces, postcodes */

  if ((status  = init_master_cfg(go, cfg, errbuf)) != eslOK)  cm_Fail(errbuf);
  if ((status  = print_run_info (go, cfg, errbuf))  != eslOK) cm_Fail(errbuf);
  
  while ((status = CMFileRead(cfg->cmfp, errbuf, &(cfg->abc), &cm)) == eslOK)
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;

      /* initialize the flags/options/params and configuration of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK)    cm_Fail(errbuf);
      print_cm_info (go, cfg, errbuf, cm);

      /* read in all sequences, this is wasteful, but Parsetrees2Alignment() requires all seqs in memory */
      seqs_to_aln = CreateSeqsToAln(100, FALSE);
      if((status = ReadSeqsToAln(cfg->abc, cfg->sqfp, 0, TRUE, seqs_to_aln, FALSE)) != eslEOF) cm_Fail("Error reading sqfile: %s\n", cfg->sqfile);
      /* align all sequences */
      if ((status = process_workunit(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);
      if ((status = output_result   (go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);
      
      /* clean up */
      FreeSeqsToAln(seqs_to_aln);
      FreeCM(cm);
      esl_sqfile_Position(cfg->sqfp, (off_t) 0); /* we may be searching this file again with another CM */
    }
  if(status != eslEOF) cm_Fail(errbuf);
  return;
}

#ifdef HAVE_MPI
/* mpi_master()
 * The MPI version of cmalign
 * Follows standard pattern for a master/worker load-balanced MPI program 
 * (SRE notes J1/78-79).
 * 
 * A master returns eslOK if it's successful. 
 * Errors in an MPI master come in two classes: recoverable and nonrecoverable.
 * If a recoverable error occurs, errbuf is filled with an error message
 * from the master or a worker, and it's sent back while returning a
 * non-eslOK error code.
 * 
 * Recoverable errors include (hopefully) all worker-side errors, and any
 * master-side error that do not affect MPI communication. Error
 * messages from recoverable messages are delayed until we've cleanly
 * shut down the workers.
 * 
 * Unrecoverable errors are master-side errors that may affect MPI
 * communication, meaning we cannot count on being able to reach the
 * workers and shut them down. Unrecoverable errors result in immediate
 * cm_Fail()'s, which will cause MPI to shut down the worker processes
 * uncleanly.
 */
static int
mpi_master(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int      xstatus       = eslOK;	/* changes from OK on recoverable error */
  int      status;
  int      have_work     = TRUE;	/* TRUE while work remains  */
  int      nproc_working = 0;	        /* number of worker processes working, up to nproc-1 */
  int      wi;          	        /* rank of next worker to get an alignment to work on */
  char    *buf           = NULL;	/* input/output buffer, for packed MPI messages */
  int      bn            = 0;
  int      pos = 1;
  int      wi_error = 0;                /* worker index that sent back an error message, if an error occurs */
  CM_t *cm;
  int nseq_per_worker;
  int nseq_this_worker;
  int nseq_prev;

  seqs_to_aln_t  *all_seqs_to_aln    = NULL;
  seqs_to_aln_t  *worker_seqs_to_aln = NULL;
  int            *seqidx         = NULL;
  long           *seedlist = NULL;       /* seeds for worker's RNGs, we send these to workers */
  
  MPI_Status mpistatus; 
  int      n;

  /* Master initialization: including, figure out the alphabet type.
   * If any failure occurs, delay printing error message until we've shut down workers.
   */
  if (xstatus == eslOK) { if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }
  if (xstatus == eslOK) { bn = 4096; if ((buf = malloc(sizeof(char) * bn))         == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seedlist       = malloc(sizeof(long) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seqidx         = malloc(sizeof(int)  * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((status         = print_run_info(go, cfg, errbuf)) != eslOK) xstatus = status; }

  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return xstatus; /* errbuf was filled above */
  ESL_DPRINTF1(("MPI master is initialized\n"));

  for (wi = 0; wi < cfg->nproc; wi++) seqidx[wi] = 0;

  for (wi = 0; wi < cfg->nproc; wi++) {
    seedlist[wi] = esl_rnd_Roll(cfg->r, 1000000000); /* not sure what to use as max for seed */
    ESL_DPRINTF1(("wi %d seed: %ld\n", wi, seedlist[wi]));
  }

  /* Worker initialization:
   * Because we've already successfully initialized the master before we start
   * initializing the workers, we don't expect worker initialization to fail;
   * so we just receive a quick OK/error code reply from each worker to be sure,
   * and don't worry about an informative message. 
   */
  for (wi = 1; wi < cfg->nproc; wi++)
    MPI_Send(&(seedlist[wi]), 1, MPI_LONG, wi, 0, MPI_COMM_WORLD);
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (status != eslOK) cm_Fail("One or more MPI worker processes failed to initialize.");
  ESL_DPRINTF1(("%d workers are initialized\n", cfg->nproc-1));
  free(seedlist);

  /* Main loop: combining load workers, send/receive, clear workers loops;
   * also, catch error states and die later, after clean shutdown of workers.
   * 
   * When a recoverable error occurs, have_work = FALSE, xstatus !=
   * eslOK, and errbuf is set to an informative message. No more
   * errbuf's can be received after the first one. We wait for all the
   * workers to clear their work units, then send them shutdown signals,
   * then finally print our errbuf and exit.
   * 
   * Unrecoverable errors just crash us out with cm_Fail().
   */

  while (xstatus == eslOK && ((status = CMFileRead(cfg->cmfp, errbuf, &(cfg->abc), &cm)) == eslOK))
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;  
      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));
      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      print_cm_info (go, cfg, errbuf, cm);
      determine_nseq_per_worker(go, cfg, cm, &nseq_per_worker); /* this func dies internally if there's some error */
      ESL_DPRINTF1(("nseq_per_worker: %d\n", nseq_per_worker));

      wi = 1;
      all_seqs_to_aln = CreateSeqsToAln(100, TRUE);
      while (have_work || nproc_working)
	{
	  if (have_work) 
	    {
	      nseq_prev = all_seqs_to_aln->nseq;
	      if((status = ReadSeqsToAln(cfg->abc, cfg->sqfp, nseq_per_worker, FALSE, all_seqs_to_aln, TRUE)) == eslOK)
		{
		  nseq_this_worker = all_seqs_to_aln->nseq - nseq_prev;
		  ESL_DPRINTF1(("MPI master read %d seqs\n", all_seqs_to_aln->nseq));
		}
	      else 
		{
		  have_work = FALSE;
		  if (status != eslEOF) cm_Fail("Sequence file read unexpectedly failed with code %d\n", status); 
		  ESL_DPRINTF1(("MPI master has run out of sequences to read (having read %d)\n", 0));
		} 
	    }
	
	  if ((have_work && nproc_working == cfg->nproc-1) || (!have_work && nproc_working > 0))
	    {
	      if (MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mpistatus) != 0) cm_Fail("mpi probe failed");
	      if (MPI_Get_count(&mpistatus, MPI_PACKED, &n)                != 0) cm_Fail("mpi get count failed");
	      wi = mpistatus.MPI_SOURCE;
	      ESL_DPRINTF1(("MPI master sees a result of %d bytes from worker %d\n", n, wi));
	      
	      if (n > bn) {
		if ((buf = realloc(buf, sizeof(char) * n)) == NULL) cm_Fail("reallocation failed");
		bn = n; 
	      }
	      if (MPI_Recv(buf, bn, MPI_PACKED, wi, 0, MPI_COMM_WORLD, &mpistatus) != 0) cm_Fail("mpi recv failed");
	      
	      /* If we're in a recoverable error state, we're only clearing worker results;
	       * just receive them, don't unpack them or print them.
	       * But if our xstatus is OK, go ahead and process the result buffer.
	       */
	      if (xstatus == eslOK) /* worker reported success. Get the result. */
		{
		  pos = 0;
		  if (MPI_Unpack(buf, bn, &pos, &xstatus, 1, MPI_INT, MPI_COMM_WORLD)     != 0)     cm_Fail("mpi unpack failed");
		  if (xstatus == eslOK) /* worker reported success. Get the results. */
		    {
		      ESL_DPRINTF1(("MPI master sees that the result buffer contains aligned sequences (seqidx: %d)\n", seqidx[wi]));
		      if ((status = cm_seqs_to_aln_MPIUnpack(cfg->abc, buf, bn, &pos, MPI_COMM_WORLD, &worker_seqs_to_aln)) != eslOK) cm_Fail("search results unpack failed");
		      ESL_DPRINTF1(("MPI master has unpacked search results\n"));
		      if ((status = add_worker_seqs_to_master(all_seqs_to_aln, worker_seqs_to_aln, seqidx[wi])) != eslOK) cm_Fail("adding worker results to master results failed");
		      /* careful not to free data from worker_seqs_to_aln we've
		       * just added to all_seqs_to_aln. we didn't copy it, we just
		       * had pointers in all_seqs_to_aln point to it. We can 
		       * free the worker's pointers to those pointers though */
		      if(worker_seqs_to_aln->sq       != NULL) free(worker_seqs_to_aln->sq);
		      if(worker_seqs_to_aln->tr       != NULL) free(worker_seqs_to_aln->tr);
		      if(worker_seqs_to_aln->cp9_tr   != NULL) free(worker_seqs_to_aln->cp9_tr);
		      if(worker_seqs_to_aln->postcode1!= NULL) free(worker_seqs_to_aln->postcode1);
		      if(worker_seqs_to_aln->postcode2!= NULL) free(worker_seqs_to_aln->postcode2);
		      if(worker_seqs_to_aln->sc       != NULL) free(worker_seqs_to_aln->sc);
		      if(worker_seqs_to_aln->pp       != NULL) free(worker_seqs_to_aln->pp);
		      if(worker_seqs_to_aln->struct_sc!= NULL) free(worker_seqs_to_aln->struct_sc);
		      free(worker_seqs_to_aln);
		    }
		  else	/* worker reported an error. Get the errbuf. */
		    {
		      if (MPI_Unpack(buf, bn, &pos, errbuf, cmERRBUFSIZE, MPI_CHAR, MPI_COMM_WORLD) != 0) cm_Fail("mpi unpack of errbuf failed");
		      ESL_DPRINTF1(("MPI master sees that the result buffer contains an error message\n"));
		      have_work = FALSE;
		      wi_error  = wi;
		    }
		}
	      nproc_working--;
	    }
	  
	  if (have_work)
	    {   
	      /* send new alignment job */
	      ESL_DPRINTF1(("MPI master is sending sequence to search to worker %d\n", wi));
	      if ((status = cm_seqs_to_aln_MPISend(all_seqs_to_aln, all_seqs_to_aln->nseq - nseq_this_worker, nseq_this_worker, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI search job send failed");
	      seqidx[wi] = all_seqs_to_aln->nseq - nseq_this_worker;
	      wi++;
	      nproc_working++;
	    }
	}
      /* if we've got valid results, output them */
      if (xstatus == eslOK) { 
	if ((status = output_result(go, cfg, errbuf, cm, all_seqs_to_aln)) != eslOK) cm_Fail(errbuf);
      }
      ESL_DPRINTF1(("MPI master: done with this CM. Telling all workers\n"));
      /* send workers the message that we're done with this CM */
      for (wi = 1; wi < cfg->nproc; wi++) 
	if ((status = cm_seqs_to_aln_MPISend(NULL, 0, 0, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("Shutting down a worker failed.");
      FreeCM(cm);
      esl_sqfile_Position(cfg->sqfp, (off_t) 0); /* we may be aligning this file again with another CM */
    }

  /* On success or recoverable errors:
   * Shut down workers cleanly. 
   */
  ESL_DPRINTF1(("MPI master is done. Shutting down all the workers cleanly\n"));
  if((cm_master_MPIBcast(NULL, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
  free(buf);
  
  if     (xstatus != eslOK) { fprintf(stderr, "Worker: %d had a problem.\n", wi_error); return xstatus; }
  else if(status != eslEOF) return status;
  else                      return eslOK; 

}

static int
mpi_worker(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int           xstatus = eslOK;
  int           status;
  CM_t         *cm  = NULL;
  char         *wbuf = NULL;	/* packed send/recv buffer  */
  int           wn   = 0;	/* allocation size for wbuf */
  int           sz, n;		/* size of a packed message */
  int           pos;
  char          errbuf[cmERRBUFSIZE];
  seqs_to_aln_t *seqs_to_aln = NULL;
  long           seed;                  /* seed for RNG, rec'd from master */
  int           i;
  MPI_Status mpistatus;           /* MPI status... */

  /* After master initialization: master broadcasts its status.
   */
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return xstatus; /* master saw an error code; workers do an immediate normal shutdown. */
  ESL_DPRINTF1(("worker %d: sees that master has initialized\n", cfg->my_rank));
  
  /* Master now sends worker initialization information (RNG seed) 
   * Workers returns their status post-initialization.
   * Initial allocation of wbuf must be large enough to guarantee that
   * we can pack an error result into it, because after initialization,
   * errors will be returned as packed (code, errbuf) messages.
   */
  if (MPI_Recv(&seed, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XEXCEPTION(eslESYS, "mpi recv failed");
  if (xstatus == eslOK) { if((cfg->r = esl_randomness_Create(seed)) == NULL)          xstatus = eslEMEM; }
  if (xstatus == eslOK) { wn = 4096;  if ((wbuf = malloc(wn * sizeof(char))) == NULL) xstatus = eslEMEM; }
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); /* everyone sends xstatus back to master */
  if (xstatus != eslOK) {
    if (wbuf != NULL) free(wbuf);
    return xstatus; /* shutdown; we passed the error back for the master to deal with. */
  }
  ESL_DPRINTF1(("worker %d: initialized seed: %ld\n", cfg->my_rank, seed));

  /* source = 0 (master); tag = 0 */
  while ((status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm)) == eslOK)
    {
      ESL_DPRINTF1(("Worker %d succesfully received CM, num states: %d num nodes: %d\n", cfg->my_rank, cm->M, cm->nodes));
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) goto ERROR;
      
      while((status = cm_seqs_to_aln_MPIRecv(cfg->abc, 0, 0, MPI_COMM_WORLD, &wbuf, &wn, &seqs_to_aln)) == eslOK)
	{
	  ESL_DPRINTF1(("worker %d: has received alignment job, nseq: %d\n", cfg->my_rank, seqs_to_aln->nseq));
	  /* align all sequences */
	  if ((status = process_workunit(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) goto ERROR;
	  ESL_DPRINTF1(("worker %d: has gathered alignment results\n", cfg->my_rank));

	  /* free the sequences, master already has a copy so we don't send them back */
	  for(i = 0; i < seqs_to_aln->nseq; i++) esl_sq_Destroy(seqs_to_aln->sq[i]);
	  free(seqs_to_aln->sq);
	  seqs_to_aln->sq = NULL;
	  n = 0;
	  if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
	    ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
	  n += sz;
	  if (cm_seqs_to_aln_MPIPackSize(seqs_to_aln, 0, seqs_to_aln->nseq, MPI_COMM_WORLD, &sz) != eslOK) 
	    ESL_XFAIL(eslFAIL, errbuf, "cm_seqs_to_aln_MPIPackSize() call failed"); 
	  n += sz;  
	  if (n > wn) {
	    void *tmp;
	    ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
	    wn = n;
	  }
	  ESL_DPRINTF1(("worker %d: has calculated the alignment results will pack into %d bytes\n", cfg->my_rank, n));
	  status = eslOK;

	  pos = 0;
	  if (MPI_Pack(&status, 1, MPI_INT, wbuf, wn, &pos, MPI_COMM_WORLD) != 0) 
	    ESL_XFAIL(eslESYS, errbuf, "mpi pack failed.");
	  if (cm_seqs_to_aln_MPIPack(seqs_to_aln, 0, seqs_to_aln->nseq, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK) 
	    ESL_XFAIL(eslFAIL, errbuf, "cm_seqs_to_aln_MPIPack() call failed"); 
	  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	  ESL_DPRINTF1(("worker %d: has sent results to master in message of %d bytes\n", cfg->my_rank, pos));
	  
	  FreeSeqsToAln(seqs_to_aln);
	}
      if(status == eslEOD) ESL_DPRINTF1(("worker %d: has seen message to stop with this CM.\n", cfg->my_rank));
      else ESL_XFAIL(eslFAIL, errbuf, "within CM loop, unexpected status code: %d received from cm_seqs_to_aln_MPIRecv()\n", status);

      FreeCM(cm);
      cm = NULL;
    }
  if (status == eslEOD) ESL_DPRINTF1(("worker %d told CMs are done.\n", cfg->my_rank));
  else ESL_FAIL(eslFAIL, errbuf, "outside CM loop, unexpected status code: %d received from cm_seqs_to_aln_MPIRecv()\n", status);

  if (wbuf != NULL) free(wbuf);
  return eslOK;

 ERROR:

  ESL_DPRINTF1(("worker %d: fails, is sending an error message, as follows:\n%s\n", cfg->my_rank, errbuf));
  pos = 0;
  MPI_Pack(&status, 1,               MPI_INT,  wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Pack(errbuf,  cmERRBUFSIZE,    MPI_CHAR, wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);

  /* if we get here this worker failed and sent an error message, now the master knows a worker
   * failed but it has to send the message to all other workers (besides this one) to abort so they 
   * can be shut down cleanly. As currently implemented, this means we have to wait here for that 
   * signal which comes in the form of a special 'empty' work packet that tells us we're done with
   * the current CM, and then a 'empty' CM broadcast that tells us we're done with all CMs in the file.
   */
  status = cm_seqs_to_aln_MPIRecv(cfg->abc, 0, 0, MPI_COMM_WORLD, &wbuf, &wn, &seqs_to_aln);
  status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm);
  /* status after each of the above calls should be eslEOD, but if it isn't we can't really do anything 
   * about it b/c we've already sent our error message, so in that scenario the MPI will break uncleanly 
   */
  return eslFAIL; /* recoverable error, master has error message and will print it */
}
#endif /*HAVE_MPI*/

static int
output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln)
{
  int status;
  ESL_MSA *msa = NULL;
  int i, ip, imax;
  float sc, struct_sc;
  void *tmp;

  /* create a new MSA, if we didn't do --inside */
  if((esl_opt_GetBoolean(go, "--cyk") || esl_opt_GetBoolean(go, "--viterbi")) || (esl_opt_GetBoolean(go, "--optacc") || esl_opt_GetBoolean(go, "--sample")))
    {
      /* optionally include a fixed alignment provided with --withali,
       * this has already been checked to see it matches the CM structure */
      if(esl_opt_GetString(go, "--withali") != NULL)
	{
	  /* grow the seqs_to_aln object */
	  if((seqs_to_aln->nseq + cfg->withmsa->nseq) > seqs_to_aln->nalloc) 
	    GrowSeqsToAln(seqs_to_aln, seqs_to_aln->nseq + cfg->withmsa->nseq - seqs_to_aln->nalloc, FALSE);
	  if((status = include_withali(go, cfg, cm, &(seqs_to_aln->sq), &(seqs_to_aln->tr), &(seqs_to_aln->nseq), errbuf)) != eslOK)
	    ESL_FAIL(status, errbuf, "--withali alignment file %s doesn't have SS_cons annotation compatible with the CM\n", esl_opt_GetString(go, "--withali"));
	}
      
      if(esl_opt_GetBoolean(go, "--viterbi"))
	{
	  ESL_DASSERT1((seqs_to_aln->cp9_tr != NULL));
	  if((status = CP9Traces2Alignment(cm, cfg->abc_out, seqs_to_aln->sq, NULL, seqs_to_aln->nseq, seqs_to_aln->cp9_tr, 
					   (! esl_opt_GetBoolean(go, "--resonly")), esl_opt_GetBoolean(go, "--matchonly"), &msa)) != eslOK)
	    goto ERROR;
	}
      else
	{
	  assert(seqs_to_aln->tr != NULL);
	  if((status = Parsetrees2Alignment(cm, cfg->abc_out, seqs_to_aln->sq, NULL, seqs_to_aln->tr, seqs_to_aln->nseq, 
					    (! esl_opt_GetBoolean(go, "--resonly")), esl_opt_GetBoolean(go, "--matchonly"), &msa)) != eslOK)
	    goto ERROR;
	}
      if(esl_opt_GetBoolean(go, "-p")) 
	{                                                                              
	  char *apostcode1;   /* aligned posterior decode array */
	  char *apostcode2;   /* aligned posterior decode array */
	  if(seqs_to_aln->postcode1 == NULL || seqs_to_aln->postcode2 == NULL) 
	    cm_Fail("-p enabled, but DispatchAlignments() did not return post codes.\n");

	  /* if --withali: we have to allocate space for the pointers to the postal codes for the preset alignment, 
	   * even though they'll be NULL. 
	   */
	  if(cfg->withmsa != NULL) {
	    ESL_RALLOC(seqs_to_aln->postcode1,tmp, sizeof(char **) * (seqs_to_aln->nseq));
	    ESL_RALLOC(seqs_to_aln->postcode2,tmp, sizeof(char **) * (seqs_to_aln->nseq));
	  }
	  imax = seqs_to_aln->nseq - 1;
	  if(cfg->withmsa != NULL) imax -= cfg->withmsa->nseq;
	  for (i = 0; i <= imax; i++)                                                   
	    {                                                                          
	      ip = (cfg->withmsa == NULL) ? i : i + cfg->withmsa->nseq;
	      if((status = make_aligned_string(msa->aseq[ip], "-_.", msa->alen, seqs_to_aln->postcode1[i], &apostcode1)) != eslOK)
		ESL_FAIL(status, errbuf, "error creating posterior string (1)\n");
	      if((status = make_aligned_string(msa->aseq[ip], "-_.", msa->alen, seqs_to_aln->postcode2[i], &apostcode2)) != eslOK)
		ESL_FAIL(status, errbuf, "error creating posterior string (2)\n");
	      esl_msa_AppendGR(msa, "POSTX.", ip, apostcode1);
	      if(! esl_opt_GetBoolean(go, "--onepost")) esl_msa_AppendGR(msa, "POST.X", ip, apostcode2);
	      free(apostcode1);                                                         
	      free(apostcode2);                                                         
	    }
	  /* if --withali: we have to allocate space for the pointers to the postal codes for the preset alignment, 
	   * even though they'll be NULL. 
	   */
	  if(cfg->withmsa != NULL) {
	    for (i = imax+1; i < seqs_to_aln->nseq; i++) 
	      seqs_to_aln->postcode1[i] = seqs_to_aln->postcode2[i] = NULL;
	  }
	}
     
#ifdef HAVE_MPI
      /* if nec, output the scores */
      if(esl_opt_GetBoolean(go, "--mpi") && (!esl_opt_GetBoolean(go, "-q"))) { 

	char *namedashes;
	int ni;
	int namewidth = 8; /* length of 'seq name' */
	/* determine the longest name in seqs_to_aln */
	for(ni = 0; ni < seqs_to_aln->nseq; ni++) namewidth = ESL_MAX(namewidth, strlen(seqs_to_aln->sq[ni]->name));
	ESL_ALLOC(namedashes, sizeof(char) * namewidth+1);
	namedashes[namewidth] = '\0';
	for(ni = 0; ni < namewidth; ni++) namedashes[ni] = '-';
	if(seqs_to_aln->struct_sc == NULL || (! NOT_IMPOSSIBLE(seqs_to_aln->struct_sc[0]))) { 
	  if(cm->align_opts & CM_ALIGN_OPTACC) { 
	    fprintf(stdout, "#\n");
	    fprintf(stdout, "# %7s  %-*s  %5s  %8s  %8s\n", "seq idx",  namewidth, "seq name",   "len",  "bit sc",   "avg prob");
	    fprintf(stdout, "# %7s  %-*s  %5s  %8s  %8s\n",  "-------", namewidth, namedashes, "-----", "--------", "--------");
	  }
	  else { 
	    fprintf(stdout, "#\n");
	    fprintf(stdout, "# %7s  %-*s  %5s  %8s\n",  "seq idx", namewidth, "seq name",   "len",  "bit sc");
	    fprintf(stdout, "# %7s  %-*s  %5s  %8s\n",  "-------", namewidth, namedashes, "-----", "--------");
	  }
	}
	else { /* we have struct scores */
	  if(cm->align_opts & CM_ALIGN_OPTACC) { 
	    fprintf(stdout, "#\n");
	    fprintf(stdout, "# %7s  %-*s  %5s  %18s  %8s\n", "",         namewidth, "",                      "",       "    bit scores    ",   "");
	    fprintf(stdout, "# %7s  %-*s  %5s  %18s  %8s\n", "",         namewidth, "",                      "",       "------------------",   "");
	    fprintf(stdout, "# %7s  %-*s  %5s  %8s  %8s  %8s\n", "seq idx",  namewidth, "seq name",   "len", "total",    "struct",   "avg prob");
	    fprintf(stdout, "# %7s  %-*s  %5s  %8s  %8s  %8s\n",  "-------", namewidth, namedashes, "-----", "--------", "--------", "--------");
	  }
	  else { 
	    fprintf(stdout, "#\n");
	    fprintf(stdout, "# %7s  %-*s  %5s  %18s\n", "", namewidth,       "",                  "",       "    bit scores    ");
	    fprintf(stdout, "# %7s  %-*s  %5s  %18s\n", "", namewidth,       "",                  "",       "------------------");
	    fprintf(stdout, "# %7s  %-*s  %5s  %8s  %8s\n",  "seq idx", namewidth,  "seq name",  "len",  "total",   "struct");
	    fprintf(stdout, "# %7s  %-*s  %5s  %8s  %8s\n",  "-------", namewidth, namedashes, "-----", "--------", "--------");
	  }
	}
	free(namedashes);

	int imin;
	float null3_correction = 0.;
	imin = (cfg->withmsa == NULL) ? 0 : cfg->withmsa->nseq;
	for (i = imin; i < seqs_to_aln->nseq; i++) {
	  ip = (cfg->withmsa == NULL) ? i : i - cfg->withmsa->nseq;
	  if(!(esl_opt_GetBoolean(go, "--no-null3"))) { ScoreCorrectionNull3CompUnknown(cm->abc, cm->null, seqs_to_aln->sq[i]->dsq, 1, seqs_to_aln->sq[i]->n, &null3_correction); }
	  if(! NOT_IMPOSSIBLE(seqs_to_aln->struct_sc[ip])) { 
	    if(cm->align_opts & CM_ALIGN_OPTACC) fprintf(stdout, "  %7d  %-*s  %5" PRId64 "  %8.2f  %8.3f\n", (ip+1), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n, seqs_to_aln->sc[ip] - null3_correction, seqs_to_aln->pp[ip]);
	    else                                 fprintf(stdout, "  %7d  %-*s  %5" PRId64 "  %8.2f\n",        (ip+1), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n, seqs_to_aln->sc[ip] - null3_correction);
	  }
	  else { /* we have struct scores */
	    if(!(esl_opt_GetBoolean(go, "--no-null3"))) seqs_to_aln->struct_sc[ip] -= ((float) ParsetreeCountMPEmissions(cm, seqs_to_aln->tr[i]) / (float) seqs_to_aln->sq[i]->n) * null3_correction; /* adjust struct_sc for NULL3 correction, this is inexact */
	    if(cm->align_opts & CM_ALIGN_OPTACC) fprintf(stdout, "  %7d  %-*s  %5" PRId64 "  %8.2f  %8.2f  %8.3f\n", (ip+1), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n, seqs_to_aln->sc[ip] - null3_correction, seqs_to_aln->struct_sc[ip], seqs_to_aln->pp[ip]);
	    else                                 fprintf(stdout, "  %7d  %-*s  %5" PRId64 "  %8.2f  %8.2f\n",        (ip+1), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n, seqs_to_aln->sc[ip] - null3_correction, seqs_to_aln->struct_sc[ip]);
	  }
	}
      }      
#endif
      if(! esl_opt_GetBoolean(go, "-q")) printf("\n");
      /* if nec, replace msa->ss_cons with ss_cons from withmsa alignment */
      if(esl_opt_GetBoolean(go, "--withpknots")) {
	if((status = add_withali_pknots(go, cfg, errbuf, cm, msa)) != eslOK) return status;
      }
      status = esl_msa_Write(cfg->ofp, msa, (esl_opt_GetBoolean(go, "-1") ? eslMSAFILE_PFAM : eslMSAFILE_STOCKHOLM));
      if      (status == eslEMEM) ESL_FAIL(status, errbuf, "Memory error when outputting alignment\n");
      else if (status != eslOK)   ESL_FAIL(status, errbuf, "Writing alignment file failed with error %d\n", status);
      
      /* if nec, output the traces */
      if(cfg->tracefp != NULL)
	{
	  for (i = 0; i < msa->nseq; i++) 
	    {
	      fprintf(cfg->tracefp, "> %s\n", seqs_to_aln->sq[i]->name);
	      if(esl_opt_GetBoolean(go,"--viterbi")) 
		{
		  fprintf(cfg->tracefp, "  SCORE : %.2f bits\n", CP9TraceScore(cm->cp9, seqs_to_aln->sq[i]->dsq, seqs_to_aln->cp9_tr[i]));
		  CP9PrintTrace(cfg->tracefp, seqs_to_aln->cp9_tr[i], cm->cp9, seqs_to_aln->sq[i]->dsq);
		}
	      else
		{
		  if((status = ParsetreeScore(cm, errbuf, seqs_to_aln->tr[i], seqs_to_aln->sq[i]->dsq, FALSE, &sc, &struct_sc)) != eslOK) return status;
		  fprintf(cfg->tracefp, "  %16s %.2f bits\n", "SCORE:", sc);
		  fprintf(cfg->tracefp, "  %16s %.2f bits\n", "STRUCTURE SCORE:", struct_sc);
		  ParsetreeDump(cfg->tracefp, seqs_to_aln->tr[i], cm, seqs_to_aln->sq[i]->dsq, NULL, NULL); /* NULLs are dmin, dmax */
		}
	      fprintf(cfg->tracefp, "//\n");
	    }
	}
      if (cfg->regressfp != NULL) {
	/* Must delete author info from msa, because it contains version
	 * and won't diff clean in regression tests. */
	if(msa->au != NULL) free(msa->au); msa->au = NULL;
	status = esl_msa_Write(cfg->regressfp, msa, eslMSAFILE_STOCKHOLM);
	if (status == eslEMEM)    ESL_FAIL(status, errbuf, "Memory error when outputting regression file\n");
	else if (status != eslOK) ESL_FAIL(status, errbuf, "Writing regression file failed with error %d\n", status);
      }
    }
  if(msa != NULL) esl_msa_Destroy(msa);
  return eslOK;

 ERROR:
  if(msa != NULL) esl_msa_Destroy(msa);
  return status;
}

/* An alignment work unit consists a seqs_to_aln_t object which contains sequences to align, 
 * and space for their parsetrees, or CP9 traces, and possibly postal codes.
 * The job is to align the sequences and create parsetrees or cp9 traces and maybe postal codes.
 */
static int
process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, 
		 seqs_to_aln_t *seqs_to_aln)
{
  int status;
  int be_quiet = esl_opt_GetBoolean(go, "-q");

#ifdef HAVE_MPI
  if(esl_opt_GetBoolean(go, "--mpi")) be_quiet = TRUE;
#endif

  if((status = DispatchAlignments(cm, errbuf, seqs_to_aln,
				  NULL, NULL, 0,  /* we're not aligning search hits */
				  esl_opt_GetInteger(go, "--banddump"),
				  esl_opt_GetInteger(go, "--dlev"), be_quiet, 
				  (! esl_opt_GetBoolean(go, "--no-null3")), cfg->r,
				  esl_opt_GetReal(go, "--mxsize"), stdout)) != eslOK) goto ERROR;
  return eslOK;
  
 ERROR:
  ESL_DPRINTF1(("worker %d: has caught an error in process_search_workunit\n", cfg->my_rank));
  FreeCM(cm);
  return status;
}

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults.
 * Configures the CM with a ConfigCM() call at end.
 */
static int
initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int status;
  int nstarts, nexits, nd;

  /* set up params/flags/options of the CM */
  cm->beta_qdb = esl_opt_GetReal(go, "--beta");
  cm->tau      = esl_opt_GetReal(go, "--tau");  /* this will be DEFAULT_TAU unless changed at command line */

  /* update cm->config_opts */
  if(esl_opt_GetBoolean(go, "-l"))
    {
      cm->config_opts |= CM_CONFIG_LOCAL;
      cm->config_opts |= CM_CONFIG_HMMLOCAL;
      cm->config_opts |= CM_CONFIG_HMMEL;
    }

  /* update cm->align_opts */
  /* optimal accuracy alignment is default */
  if(esl_opt_GetBoolean(go, "--optacc"))      cm->align_opts  |= CM_ALIGN_OPTACC;
  if(esl_opt_GetBoolean(go, "--sample"))      cm->align_opts  |= CM_ALIGN_SAMPLE;
  if(esl_opt_GetBoolean(go, "--hbanded"))     cm->align_opts  |= CM_ALIGN_HBANDED;
  if(esl_opt_GetBoolean(go, "--nonbanded"))   cm->align_opts  &= ~CM_ALIGN_HBANDED;
  if(esl_opt_GetBoolean(go, "--sub"))         cm->align_opts  |= CM_ALIGN_SUB;
  if(esl_opt_GetBoolean(go, "--viterbi"))     cm->align_opts  |= CM_ALIGN_HMMVITERBI;
  if(esl_opt_GetBoolean(go, "--small"))       cm->align_opts  |= CM_ALIGN_SMALL;
  if(esl_opt_GetBoolean(go, "-p"))            cm->align_opts  |= CM_ALIGN_POST;
  if(esl_opt_GetBoolean(go, "--hsafe"))       cm->align_opts  |= CM_ALIGN_HMMSAFE;
  if(esl_opt_GetBoolean(go, "--fins"))        cm->align_opts  |= CM_ALIGN_FLUSHINSERTS;

  if(esl_opt_GetBoolean(go, "--inside"))      cm->align_opts  |= CM_ALIGN_INSIDE;
  if(esl_opt_GetBoolean(go, "--checkpost"))   cm->align_opts  |= CM_ALIGN_CHECKINOUT;
  if(esl_opt_GetBoolean(go, "--checkfb"))     cm->align_opts  |= CM_ALIGN_CHECKFB;
  if(esl_opt_GetBoolean(go, "--sums"))        cm->align_opts  |= CM_ALIGN_SUMS;
  /* config QDB? */
  if(esl_opt_GetBoolean(go, "--qdb"))          
    { 
      cm->align_opts  |= CM_ALIGN_QDB;
      cm->config_opts |= CM_CONFIG_QDB;
    }

  /* BEGIN (POTENTIALLY) TEMPORARY BLOCK */
  /* set aggregate local begin/end probs, set with --pbegin, --pend, defaults are DEFAULT_PBEGIN, DEFAULT_PEND */
  cm->pbegin = esl_opt_GetReal(go, "--pbegin");
  cm->pend   = esl_opt_GetReal(go, "--pend");
  /* possibly overwrite local begin probs such that all begin points are equiprobable (--pebegin) */
  if(esl_opt_GetBoolean(go, "--pebegin")) {
    nstarts = 0;
    for (nd = 2; nd < cm->nodes; nd++) 
      if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd || cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd) 
	nstarts++;
    /* printf("nstarts: %d\n", nstarts); */
    cm->pbegin = 1.- (1./(1+nstarts));
    /* printf("pbegin: %.5f\n", cm->pbegin); */
  }
  /* possibly overwrite cm->pend so that local end prob from all legal states is fixed,
   * this is strange in that cm->pend may be placed as a number greater than 1., this number
   * is then divided by nexits in ConfigLocalEnds() to get the prob for each v --> EL transition,
   * this is guaranteed by the way we calculate it to be < 1.,  it's the argument from --pfend */
  if(! esl_opt_IsDefault(go, "--pfend")) {
    nexits = 0;
    for (nd = 1; nd < cm->nodes; nd++) {
      if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	   cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	   cm->ndtype[nd] == BEGR_nd) && 
	  cm->ndtype[nd+1] != END_nd)
	nexits++;
    }
    cm->pend = nexits * esl_opt_GetReal(go, "--pfend");
  }
  /* END (POTENTIALLY) TEMPORARY BLOCK */


  /* finally, configure the CM for alignment based on cm->config_opts and cm->align_opts.
   * set local mode, make cp9 HMM, calculate QD bands etc. 
   */
  ConfigCM(cm, FALSE);  /* FALSE says do not calculate W unless nec b/c we're using QDBs */

  /* if(cfg->my_rank == 0) printf("CM %d: %s\n", (cfg->ncm), cm->name); 
   * debug_print_cm_params(stdout, cm);
   */

  /* if we're master and we're trying to include an alignment, make sure it is consistent with CM structure */
  if(cfg->my_rank == 0 && (! esl_opt_IsDefault(go, "--withali"))) { 
    if((status = check_withali(go, cfg, cm, &(cfg->withmsa), &(cfg->withali_mtr))) != eslOK)
      ESL_FAIL(status, errbuf, "--withali alignment file %s doesn't have a SS_cons compatible with the CM\n", esl_opt_GetString(go, "--withali"));
    }
  return eslOK;
}

/* Function: compare_cm_guide_trees()
 * EPN, Tue Mar  6 08:32:12 2007
 *
 * Purpose:  Given two CMs, cm1 and cm2, compare them, returning TRUE 
 *           iff they have the same guide tree (same node architecture).
 *
 * Args:     cm1          - covariance model number 1
 *           cm2          - covariance model number 2
 * 
 * Returns:  TRUE if CMs have same guide tree, FALSE otherwise
 */
static int compare_cm_guide_trees(CM_t *cm1, CM_t *cm2)
{
  int          nd; 
  if(cm1->nodes != cm2->nodes) return FALSE;
  for(nd = 0; nd < cm1->nodes; nd++)
    if(cm1->ndtype[nd] != cm2->ndtype[nd]) return FALSE;
  return TRUE;
}

/* Function: check_withali()
 * EPN, Tue Mar  6 06:25:02 2007
 *
 * Purpose:  Ensure that the alignment to include has a secondary
 *           structure that matches our CM. Pass the alignment back
 *           as *ret_msa.
 *
 * Returns:  <eslOK> on success.
 *           <eslEINCOMPAT> if alignment doesn't match the CM 
 */
static int check_withali(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, ESL_MSA **ret_msa, Parsetree_t **ret_mtr)
{
  int           status;
  ESL_MSA      *msa      = NULL; /* alignment we're including  */
  CM_t         *new_cm   = NULL; /* CM built from MSA, we check it has same guide tree as 'cm' */
  Parsetree_t  *mtr      = NULL; /* master structure tree from the alignment*/
  char          errbuf[cmERRBUFSIZE];

  /* cfg->withalifp is open */
  status = esl_msa_Read(cfg->withalifp, &msa);
  if (status == eslEFORMAT)  cm_Fail("--withali alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", cfg->withalifp->linenumber, cfg->withalifp->fname, cfg->withalifp->errbuf, cfg->withalifp->buf);
  else if (status != eslOK)       cm_Fail("--withali alignment file read unexpectedly failed with code %d\n", status);

  /* Some input data cleaning. */
  if (esl_opt_GetBoolean(go, "--rf") && msa->rf == NULL) 
    ESL_FAIL(eslFAIL, errbuf, "--rf invoked but --withali alignment has no reference coord annotation.\n");
  if (msa->ss_cons == NULL) 
    ESL_FAIL(eslFAIL, errbuf, "--withali alignment did not contain consensus structure annotation.\n");
  if (esl_opt_GetBoolean(go, "--withpknots")) /* copy the original secondary structure */
    esl_strdup(msa->ss_cons, -1, &(cfg->withss_cons));
 if (! clean_cs(msa->ss_cons, msa->alen, TRUE))
    ESL_FAIL(eslFAIL, errbuf, "Failed to parse consensus structure annotation for --withali alignment\n");

  /* Build a CM from a master guide tree built from the msa, 
   * then check to make sure this CM has same emit map as the CM
   * we've had passed in. This is fragile and hopefully temporary. 
   * Another solution would be to use a checksum, but CM files don't 
   * have checksums yet.
   */
  if((status = HandModelmaker(msa, errbuf, esl_opt_GetBoolean(go, "--rf"), esl_opt_GetReal(go, "--gapthresh"), &new_cm, &mtr)) != eslOK) return status;
  if(!(compare_cm_guide_trees(cm, new_cm))) {
    /* no need to try rebalancing, that doesn't change the guidetree (seriously this is from cm.c::CMRebalance():
     * for (x = 0; x < new->nodes; x++) new->ndtype[x]  = cm->ndtype[x];
     */
    status = eslEINCOMPAT;
    goto ERROR;
  }

  /* if we get here, the CM guide trees match */
  if(new_cm   != NULL) FreeCM(new_cm);
  *ret_mtr = mtr;
  *ret_msa = msa;
  return eslOK;

 ERROR:
  if(msa != NULL)      esl_msa_Destroy(msa);
  if(new_cm   != NULL) FreeCM(new_cm);
  if(mtr != NULL)      FreeParsetree(mtr);
  return eslEINCOMPAT;
}

/* Function: include_withali()
 * EPN, Tue Mar  6 06:25:02 2007
 *
 * Purpose:  Determine the implicit parses of sequences in an
 *           MSA to a CM and append them to passed in data structures.
 *           We've already checked to make sure the MSA's consensus
 *           structure matches the CM.
 *
 * Args:     go           - command line options
 *           cfg          - cmalign configuration, includes msa to add
 *           cm           - CM we're aligning to 
 *           ret_sq       - pre-existing sequences, to append msa seqs to
 *           ret_tr       - pre-existing parsetrees, to append to
 *           ret_nseq     - number of exisiting seqs, updated at end of function
 *           errbuf       - easel error message
 * 
 * Returns:  eslOK on success, eslEMEM on memory error
 *           Also new, realloc'ed arrays for sq, tr in ret_seq, ret_tr; 
 *           *ret_nseq is increased by number of seqs in cfg->withmsa.
*/
static int include_withali(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr, int *ret_nseq, char *errbuf)
{
  int           status;
  void         *tmp;      /* for ESL_RALLOC() */
  int           i;	  /* counter over aseqs       */
  int           ip;	  /* offset counter over aseqs */
  char        **uaseq;    /* unaligned seqs, dealigned from the MSA */
  char         **aseq;    /*   aligned text seqs */
  int           apos;     /*   aligned position index */
  int           uapos;    /* unaligned position index */
  int           x;        /* counter of parsetree nodes */
  int         **map;      /* [0..msa->nseq-1][0..msa->alen] map from aligned
			   * positions to unaligned (non-gap) positions */


  /* Contract check */
  if(cfg->withmsa == NULL)                     cm_Fail("ERROR in include_withali() withmsa is NULL.\n");
  if(! (cfg->withmsa->flags & eslMSA_DIGITAL)) cm_Fail("ERROR in include_withali() withmsa is not digitized.\n");

  /* For each seq in the MSA, map the aligned sequences coords to 
   * the unaligned coords, we stay in digitized seq coords (1..alen),
   * we need this for converting parsetrees from Transmogrify (which
   * have emitl and emitr in aligned coords) to unaligned coords, so 
   * we can call Parsetrees2Alignment() with them. */
  ESL_ALLOC(map,   sizeof(int *)  * cfg->withmsa->nseq);
  ESL_ALLOC(uaseq, sizeof(char *) * cfg->withmsa->nseq);
  ESL_ALLOC(aseq,  sizeof(char *) * cfg->withmsa->nseq);
  for (i = 0; i < cfg->withmsa->nseq; i++)
    {
      ESL_ALLOC(map[i],   sizeof(int)  * (cfg->withmsa->alen+1));
      ESL_ALLOC(aseq[i],  sizeof(char) * (cfg->withmsa->alen+1));
      map[i][0] = -1; /* invalid */
      uapos = 1;
      for(apos = 1; apos <= cfg->withmsa->alen; apos++)
	{
	  if (!esl_abc_XIsGap(cfg->withmsa->abc, cfg->withmsa->ax[i][apos]))
	    map[i][apos] = uapos++;
	  else
	    map[i][apos] = -1;
	}
      /* we need digitized AND text seqs for Transmogrify */
      esl_abc_Textize(cfg->withmsa->abc, cfg->withmsa->ax[i], cfg->withmsa->alen, aseq[i]);
      esl_strdup(aseq[i], -1, &(uaseq[i]));
      esl_strdealign(uaseq[i], uaseq[i], "-_.", NULL);
    }
  ESL_RALLOC((*ret_tr),  tmp, (sizeof(Parsetree_t *)  * (*ret_nseq + cfg->withmsa->nseq)));
  ESL_RALLOC((*ret_sq),  tmp, (sizeof(ESL_SQ *)       * (*ret_nseq + cfg->withmsa->nseq)));

  /* Transmogrify each aligned seq to get a parsetree */
  /*for (i = 0; i < cfg->withmsa->nseq; i++)*/
  for (i = *ret_nseq; i < (*ret_nseq + cfg->withmsa->nseq); i++)
    {
      ip = i - *ret_nseq;
      (*ret_tr)[i] = Transmogrify(cm, cfg->withali_mtr, cfg->withmsa->ax[ip], aseq[ip], cfg->withmsa->alen);
      /* ret_tr[i] is in alignment coords, convert it to unaligned coords, */
      for(x = 0; x < (*ret_tr)[i]->n; x++)
	{
	  if((*ret_tr)[i]->emitl[x] != -1)
	    (*ret_tr)[i]->emitl[x] = map[ip][(*ret_tr)[i]->emitl[x]];
	  if((*ret_tr)[i]->emitr[x] != -1)
	    (*ret_tr)[i]->emitr[x] = map[ip][(*ret_tr)[i]->emitr[x]];
	}
      (*ret_sq)[i]      = esl_sq_CreateFrom(cfg->withmsa->sqname[ip], uaseq[ip], NULL, NULL, NULL);
      esl_sq_Digitize(cm->abc, (*ret_sq)[i]);
    }

  /* Swap some pointers so the included alignment appears at the top of the output 
   * alignment instead of the bottom. */
  Parsetree_t **tmp_tr;
  ESL_SQ      **tmp_sq;
  ESL_ALLOC(tmp_tr, sizeof(Parsetree_t *) * (*ret_nseq + cfg->withmsa->nseq));
  ESL_ALLOC(tmp_sq, sizeof(ESL_SQ *)      * (*ret_nseq + cfg->withmsa->nseq));
  for(i = 0; i < (*ret_nseq + cfg->withmsa->nseq); i++)
    {
      tmp_tr[i] = (*ret_tr)[i];
      tmp_sq[i] = (*ret_sq)[i];
    }
  for(i = 0; i < *ret_nseq; i++)
    {
      ip = i + cfg->withmsa->nseq;
      (*ret_tr)[ip] = tmp_tr[i];
      (*ret_sq)[ip] = tmp_sq[i];
    }
  for(i = *ret_nseq; i < (*ret_nseq + cfg->withmsa->nseq); i++)
    {
      ip = i - *ret_nseq;
      (*ret_tr)[ip] = tmp_tr[i];
      (*ret_sq)[ip] = tmp_sq[i];
    }
  free(tmp_tr);
  free(tmp_sq);

  /* update *ret_nseq */
  *ret_nseq    += cfg->withmsa->nseq;

  /* Clean up and exit. */
  esl_Free2D((void **) map,   cfg->withmsa->nseq);
  esl_Free2D((void **) uaseq, cfg->withmsa->nseq);
  esl_Free2D((void **) aseq,  cfg->withmsa->nseq);
  return eslOK;

 ERROR:
  esl_Free2D((void **) map,   cfg->withmsa->nseq);
  esl_Free2D((void **) uaseq, cfg->withmsa->nseq);
  return status;
}

/* Function: make_aligned_string() 
 * Incept:   EPN, Thu Aug  2 18:24:49 2007
 *           stolen from Squid:alignio.c:MakeAlignedString()
 * 
 * Purpose:  Given a raw string of some type (secondary structure, say),
 *           align it to a given aseq by putting gaps wherever the
 *           aseq has gaps. 
 *           
 * Args:     aseq:  template for alignment
 *           gapstring: defines all gap chars ex: "-_."
 *           alen:  length of aseq
 *           ss:    raw string to align to aseq
 *           ret_s: RETURN: aligned ss
 *           
 * Return:   eslOK on success, eslEMEM on memory allocation error,
 *           eslEINCONCEIVABLE on strange error,
 *           ret_ss is alloc'ed here and must be free'd by caller.
 */
int
make_aligned_string(char *aseq, char *gapstring, int alen, char *ss, char **ret_s)
{
  int status;
  char *new; 
  int   apos, rpos;
  int   rlen;

  ESL_ALLOC(new, (sizeof(char) * (alen+1)));
  for (apos = rpos = 0; apos < alen; apos++)
    {
      if (strchr(gapstring, aseq[apos]) != NULL)
	new[apos] = aseq[apos];
      else
	new[apos] = ss[rpos++];
    }
  new[apos] = '\0';
  
  rlen = strlen(ss);
  if (rpos != rlen)
    {
      if(new != NULL) free(new);
      return eslEINCONCEIVABLE;
    }
  *ret_s = new;
  return eslOK;

 ERROR:
  if(new != NULL) free(new);
  return status;
}

/* Function: add_withali_pknots()
 * EPN, Wed Oct 17 18:24:33 2007
 *
 * Purpose:  Determine the pseudoknots that were in consensus columns of
 *           the --withali alignment and add them to newmsa->ss_cons.
 *
 * Args:     go           - command line options
 *           cfg          - cmalign configuration, includes msa to add
 *           errbuf       - easel error message
 *           cm           - the CM, only used to get cm->clen 
 *           newmsa       - MSA from Parsetrees2Alignment(), we want to add to it's ss_cons
 * 
 * Returns:  eslOK on success, eslEMEM on memory error
 *           eslEINVAL on unpredicted situation
*/
static int add_withali_pknots(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_MSA *newmsa)
{
  int           status;
  int           apos;     /*   aligned position index */
  int           cpos;     /* consensus column index */
  int           ngaps;
  float         gapthresh;
  int           withmsa_clen = 0;
  int           idx;     /* sequence index */
  int           i, j, i_cpos, j_cpos; /* residue position indices */
  /* Contract check */
  if(cfg->withmsa == NULL) cm_Fail("ERROR in add_withali_pknots() cfg->withmsa is NULL.\n");
  if(cfg->withss_cons  == NULL) cm_Fail("ERROR in add_withali_pknots() cfg->withss_cons is NULL.\n");
  if(! (cfg->withmsa->flags & eslMSA_DIGITAL)) cm_Fail("ERROR in add_withali_pknots() cfg->withmsa is not digitized.\n");

  /* 10 easy, convoluted steps. One reason for so many steps is 
   * we can't build ss_cons strings from pseudoknotted ct arrays,
   * so we have to work around that. 
   */

  /* 1. determine consensus columns of withmsa.  
   * If we've gotten this far, there should be same number as
   * cm->clen. code block stolen from cm_modelmaker.c matassign is
   * 1..alen. Values are 1 if a match (consensus) column, 0 if insert
   * column.
   */
  gapthresh = esl_opt_GetReal(go, "--gapthresh");
  int *matassign;
  ESL_ALLOC(matassign, sizeof(int) * (cfg->withmsa->alen+1));
  /* Watch for off-by-one. rf is [0..alen-1]; matassign is [1..alen] */
  if (esl_opt_GetBoolean(go, "--rf")) {
    for (apos = 1; apos <= cfg->withmsa->alen; apos++) 
      matassign[apos] = (esl_abc_CIsGap(cfg->withmsa->abc, cfg->withmsa->rf[apos-1])? FALSE : TRUE);
  }
  else { /* --rf not enabled */
    for (apos = 1; apos <= cfg->withmsa->alen; apos++) {
      for (ngaps = 0, idx = 0; idx < cfg->withmsa->nseq; idx++)
	if (esl_abc_XIsGap(cfg->withmsa->abc, cfg->withmsa->ax[idx][apos])) ngaps++;
      matassign[apos] = ((double) ngaps / (double) cfg->withmsa->nseq > gapthresh) ? 0 : 1;
    }
  }
  for (apos = 1; apos <= cfg->withmsa->alen; apos++) withmsa_clen += matassign[apos];
  if(withmsa_clen != cm->clen) ESL_FAIL(eslFAIL, errbuf, "withmsa consensus length != cm consensus length. A previous check for this passed, this is a coding error.");

  /* 2. get ct array for consensus structure of withmsa BEFORE we stripped away it's pknots,
   * this was saved in check_withali().
   */
  int *ct;
  ESL_ALLOC(ct, (cfg->withmsa->alen+1) * sizeof(int));
  if (esl_wuss2ct(cfg->withss_cons, cfg->withmsa->alen, ct) != eslOK)  
    ESL_FAIL(eslFAIL, errbuf, "withmsa original ss_cons inconsistent. A previous check for this passed, this is a coding error.");

  /* 3. also get a ct with no pknots, we'll need this to figure out where the pknots go */
  int *ct_noknots;
  ESL_ALLOC(ct_noknots, (cfg->withmsa->alen+1) * sizeof(int));
  if (esl_wuss2ct(cfg->withmsa->ss_cons, cfg->withmsa->alen, ct_noknots) != eslOK)  
    ESL_FAIL(eslFAIL, errbuf, "withmsa original ss_cons inconsistent. A previous check for this passed, this is a coding error.");

  /* 4. remove any base pairs (i,j) from ct and ct_noknots for which i or j are non-consensus columns */
  for (apos = 1; apos <= cfg->withmsa->alen; apos++) {
    if(! matassign[apos]) { /* apos is not a consensus column */
      if(ct[apos] != 0) ct[ct[apos]] = 0;
      ct[apos] = 0;
      if(ct_noknots[apos] != 0) ct_noknots[ct_noknots[apos]] = 0;
      ct_noknots[apos] = 0;
    }
  }

  /* 5. get a map from alignment coords to consensus coords */
  int *a2c_map;
  ESL_ALLOC(a2c_map, sizeof(int) * (cfg->withmsa->alen + 1));
  a2c_map[0] = -1;
  cpos = 1;
  for(apos = 1; apos <= cfg->withmsa->alen; apos++) {
    if(matassign[apos]) a2c_map[apos] = cpos++;
    else a2c_map[apos] = 0;
  }

  /* 6. use a2c_map to create c_ct_noknots, which is just ct_noknots
   * changed from aligned coordinates to consensus column
   * coordinates. remember no non-consensus column cpos should have
   * ct_noknots[cpos] != 0, because we stripped those bps.
   */
  int *c_ct_noknots;
  ESL_ALLOC(c_ct_noknots, sizeof(int) * (cm->clen+1));
  esl_vec_ISet(c_ct_noknots, (cm->clen+1), 0);
  for(apos = 1; apos <= cfg->withmsa->alen; apos++) {
    i = apos; j = ct_noknots[i];
    if(j != 0 && i < j) { /* if i > j, we've already covered it */
      if(a2c_map[i] == 0) ESL_FAIL(eslFAIL, errbuf, "withmsa apos: %d has structure, but is not consensus, this should never happen.", i);
      if(a2c_map[j] == 0) ESL_FAIL(eslFAIL, errbuf, "withmsa apos: %d has structure, but is not consensus, this should never happen.", j);
      c_ct_noknots[a2c_map[i]] = a2c_map[j];
      c_ct_noknots[a2c_map[j]] = a2c_map[i];
    }      
  }
  for(cpos = 1; cpos <= cm->clen; cpos++)
    printf("ct[%d]: %d\n", cpos, c_ct_noknots[cpos]);
  
  /* 7. build new consensus structure, with no knots and only consensus columns */
  char *c_sscons;
  ESL_ALLOC(c_sscons, sizeof(char) * (cm->clen + 1));
  if((status = esl_ct2wuss(c_ct_noknots, cm->clen, c_sscons)) != eslOK) cm_Fail("ct2wuss failed with (supposedly) no knots");
  
  /* 8. add back in consensus knots, using a2c_map */
  for(apos = 1; apos <= cfg->withmsa->alen; apos++) {
    if(matassign[apos]) {
      i = apos; j = ct[i];
      if(ct[i] != 0 && i < j) { /* if i > j, we've already updated it */
	i_cpos = a2c_map[i];
	j_cpos = a2c_map[j];
	/* printf("\t\tct bp i: %3d i_cpos: %3d j_cpos: %d j_cpos: %3d\n", i, i_cpos, j, j_cpos);
	   printf("\t\tc_sscons[i_cpos-1]: %c\n", c_sscons[(i_cpos-1)]);
	   printf("\t\tc_sscons[j_cpos-1]: %c\n", c_sscons[(j_cpos-1)]);
	   printf("\t\tcfg->withss_cons[(i-1)]: %c\n", cfg->withss_cons[(i-1)]);
	   printf("\t\tcfg->withss_cons[(j-1)]: %c\n", cfg->withss_cons[(j-1)]);
	*/
	c_sscons[(i_cpos-1)] = cfg->withss_cons[(i-1)]; /* add pknot annotation for left bp */
	c_sscons[(j_cpos-1)] = cfg->withss_cons[(j-1)]; /* add pknot annotation for right bp */
      }
    }
  }

  /* 9. c_sscons is the pknotted structure, but limited to the consensus columns,
   * final step is to overwrite newmsa->ss_cons characters for consensus columns 
   * only, by simply replacing them with characters from c_sscons.
   * we need a new map from consensus columns to newmsa align coords first,
   * we use the fact that newmsa->rf columns that are non-gapped are consensus
   * columns.
   */
  int *new_c2a_map;
  ESL_ALLOC(new_c2a_map, sizeof(int) * (cm->clen + 1));
  new_c2a_map[0] = -1;
  new_c2a_map[0] = -1;
  cpos = 0;
  for(apos = 1; apos <= newmsa->alen; apos++) 
    if(! esl_abc_CIsGap(newmsa->abc, newmsa->rf[(apos-1)])) new_c2a_map[++cpos] = apos;

  /* 10. overwrite newmsa->ss_cons */
  for(cpos = 1; cpos <= cm->clen; cpos++) 
    newmsa->ss_cons[(new_c2a_map[cpos]-1)] = c_sscons[(cpos-1)];
  
  /* free memory and return */
  free(new_c2a_map);
  free(c_sscons);
  free(a2c_map);
  free(ct_noknots);
  free(c_ct_noknots);
  free(ct);
  free(matassign);
  return eslOK;

 ERROR:
  return status;
}

/* Function: print_run_info
 * Date:     EPN, Thu Feb 28 14:26:42 2008
 *
 * Purpose:  Print information on this run of cmalign.
 *           Command used to run it, and execution date.
 *
 * Returns:  eslOK on success
 */
static int
print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf)
{
  int status;
  char *command;
  char *date;

  if(esl_opt_GetBoolean(go, "-q")) return eslOK;

  if((status = get_command(go, errbuf, &command)) != eslOK) return status;
  if((status = GetDate    (errbuf, &date))        != eslOK) return status;

  fprintf(stdout, "%-10s %s\n",  "# command:", command);
  fprintf(stdout, "%-10s %s\n",  "# date:",    date);
  if(cfg->nproc > 1) fprintf(stdout, "# %-8s %d\n", "nproc:", cfg->nproc);
  if(esl_opt_GetBoolean(go, "--sample")) fprintf(stdout, "%-10s %ld\n", "# seed:", esl_randomness_GetSeed(cfg->r));

  fprintf(stdout, "#\n");
  free(command);
  free(date);
  return eslOK;
}

/* Function: print_cm_info
 * Date:     EPN, Thu Feb 28 14:44:17 2008
 *
 * Purpose:  Print per-CM info to output file (stdout unless -o). 
 */
static void
print_cm_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm)
{

  if(esl_opt_GetBoolean(go, "-q")) return;

  int do_hbanded = (cm->align_opts & CM_ALIGN_HBANDED) ? TRUE : FALSE;
  int do_qdb     = (cm->align_opts & CM_ALIGN_QDB)     ? TRUE : FALSE;
  if(cm->align_opts & CM_ALIGN_HMMVITERBI) { do_hbanded = do_qdb = FALSE; }

  fprintf(stdout, "# %-25s  %9s  %6s  %3s  %5s  %6s\n", "cm name",                   "algorithm", "config", "sub", "bands", (do_hbanded) ? "tau" : ((do_qdb) ? "beta" : "")); 
  fprintf(stdout, "# %-25s  %9s  %6s  %3s  %5s  %6s\n", "-------------------------", "---------", "------", "---", "-----", (do_hbanded || do_qdb) ? "------" : ""); 
  fprintf(stdout, "# %-25.25s  %9s  %6s  %3s", 
	  cm->name,
	  ((esl_opt_GetBoolean(go, "--cyk")) ? "cyk" : ((esl_opt_GetBoolean(go, "--sample")) ? "sample" : (esl_opt_GetBoolean(go, "--viterbi") ? "hmm vit" : "opt acc"))), 
	  (esl_opt_GetBoolean(go, "-l")) ? "local" : "global",
	  (esl_opt_GetBoolean(go, "--sub")) ? "yes" : "no");
  /* bands and beta/tau */
  if     (do_hbanded)    fprintf(stdout, "  %5s  %6.0e\n", "hmm", cm->tau);
  else if(do_qdb)        fprintf(stdout, "  %5s  %6.0e\n", "qdb", cm->beta_qdb);
  else                   fprintf(stdout, "  %5s  %6s\n", "none", "");

  return;
}

/* Function: get_command
 * Date:     EPN, Fri Jan 25 13:56:10 2008
 *
 * Purpose:  Return the command used to call cmscore
 *           in <ret_command>.
 *
 * Returns:  eslOK on success; eslEMEM on allocation failure.
 */
int 
get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command)
{
  int status;
  int i;
  char *command = NULL;

  for (i = 0; i < go->argc; i++) { /* copy all command line options and args */
    if((status = esl_strcat(&(command),  -1, go->argv[i], -1)) != eslOK) goto ERROR;
    if(i < (go->argc-1)) if((status = esl_strcat(&(command), -1, " ", 1)) != eslOK) goto ERROR;
  }
  *ret_command = command;

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "get_command(): memory allocation error.");
  return status;
}

/* serial_merge_alignments_only()
 *
 * A special mode of cmalign invoked with --merge <f>. In this case
 * we do not align sequences to a CM, rather we only merge two alignments
 * created by previous runs of cmalign (to the same CM) together. 
 * We output the alignment to either stdout or file <of> of -o <of> 
 * and return, the program will then skip serial_master() and exit.
 * If we run into any problem here we cm_Fail() with an informative
 * error message. --merge is incompatible with --mpi, so we don't have to
 * worry about MPI mode with --merge.
 */
static void
serial_merge_alignments_only(const ESL_GETOPTS *go)
{
  int status;
  char errbuf[cmERRBUFSIZE];
  CMFILE       *cmfp;		/* open input CM file stream       */
  ESL_MSAFILE  *msafp1;         /* open alifile 1 (argv[2]) */
  ESL_MSAFILE  *msafp2;         /* open alifile 2 (from --merge <f>) */
  char *cmfile;
  char *msafile1;
  char *msafile2;
  ESL_MSA *msa1;
  ESL_MSA *msa2;
  ESL_MSA *merged_msa;
  CM_t *cm, *tmp_cm;
  ESL_ALPHABET *abc, *abc_out;
  abc = NULL;
  Parsetree_t *mtr1, *mtr2;
  FILE *ofp;
  int nmerged_tr, nmerged_sq;
  int msa1_nseq, msa2_nseq;
  int i;
  Parsetree_t **tr1, **tr2;
  ESL_SQ **sq1, **sq2;
  
  /* markup from msa1 and msa2, we'll add it to the merged_msa */
  int ncomment1, ncomment2, ngf1, ngf2, ngs1, ngs2, ngr1, ngr2;
  char **gf_tag1, **gf_tag2, **gs_tag1, **gs_tag2, **gr_tag1, **gr_tag2;
  char  **comment1, **comment2, **gf1, **gf2;
  char ***gs1, ***gs2, ***gr1, ***gr2;

  /* open output file */
  if (esl_opt_GetString(go, "-o") != NULL) {
    if ((ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
      cm_Fail("Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
  } else ofp = stdout;

  /* open CM file, read first CM */
  cmfile = esl_opt_GetArg(go, 1); 
  if((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if      (esl_opt_GetBoolean(go, "--rna")) abc_out = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna")) abc_out = esl_alphabet_Create(eslDNA);
  if((status = CMFileRead(cmfp, errbuf, &abc, &cm)) != eslOK) cm_Fail(errbuf);
  /* configure the CM, only really nec so we get the appropriate #=GC RF line in the output alignment
   * which is built in Parsetrees2Alignment() using the log-odds emissions scores
   */
  ConfigCM(cm, FALSE);  /* FALSE says do not calculate W unless nec b/c we're using QDBs */

  /* open msa file 1, argv[2] */
  msafile1 = esl_opt_GetArg(go, 2); 
  status = esl_msafile_Open(msafile1, eslMSAFILE_UNKNOWN, NULL, &msafp1);
  if      (status == eslENOTFOUND) cm_Fail("Alignment file %s doesn't exist or is not readable\n", msafile1);
  else if (status == eslEFORMAT)   cm_Fail("Couldn't determine format of alignment %s\n", msafile1);
  else if (status != eslOK)        cm_Fail("Alignment file %s open failed with error %d\n", msafile1, status);
  /* read first alignment from msa file 1 */
  esl_msafile_SetDigital(msafp1, abc);
  if((status = esl_msa_Read(msafp1, &msa1)) != eslOK) cm_Fail("No alignments in %s\n", msafile1);

  /* open msa file 2, argv[3] */
  msafile2 = esl_opt_GetArg(go, 3);
  status = esl_msafile_Open(msafile2, eslMSAFILE_UNKNOWN, NULL, &msafp2);
  if      (status == eslENOTFOUND) cm_Fail("Alignment file %s doesn't exist or is not readable\n", msafile2);
  else if (status == eslEFORMAT)   cm_Fail("Couldn't determine format of alignment %s\n", msafile2);
  else if (status != eslOK)        cm_Fail("Alignment file %s open failed with error %d\n", msafile2, status);
  /* read first alignment from msa file 2 */
  esl_msafile_SetDigital(msafp2, cm->abc);
  if((status = esl_msa_Read(msafp2, &msa2)) != eslOK) cm_Fail("No alignments in %s\n", msafile2);

  /* validation: Build CMs from both MSAs, the guidetrees should match the guidetree for the CM we read from the cmfile */
  if((status = HandModelmaker(msa1, errbuf, TRUE, 0.5, &tmp_cm, &(mtr1))) != eslOK) cm_Fail(errbuf);
  /*                                        !use RF! */
  if(!(CompareCMGuideTrees(cm, tmp_cm))) { 
    /* no need to try rebalancing, that doesn't change the guidetree (seriously this is from cm.c::CMRebalance():
     * for (x = 0; x < new->nodes; x++) new->ndtype[x]  = cm->ndtype[x];
     */
    cm_Fail("with --merge, alignment in %s could not have been created using the first CM in %s.", msafile1, cmfile);
  }
  FreeCM(tmp_cm);

  if((status = HandModelmaker(msa2, errbuf, TRUE, 0.5, &tmp_cm, &(mtr2))) != eslOK) cm_Fail(errbuf);
  /*                                        !use RF! */
  if(!(CompareCMGuideTrees(cm, tmp_cm))) { 
    /* no need to try rebalancing, that doesn't change the guidetree (seriously this is from cm.c::CMRebalance():
     * for (x = 0; x < new->nodes; x++) new->ndtype[x]  = cm->ndtype[x];
     */
    cm_Fail("with --merge, alignment in %s could not have been created using the first CM in %s.", msafile2, cmfile);
  }
  FreeCM(tmp_cm);
  /* now we know msa1 and msa2 both have SS_cons lines that match the first CM we read from cmfile, so they must match each other */

  /* store the unparsed stockholm markup from each msa, we'll regurgitate it in the new alignment */
  if((status = get_markup_from_msa(msa1, errbuf, &comment1, &ncomment1, &gf_tag1, &gf1, &ngf1, &gs_tag1, &gs1, &ngs1, &gr_tag1, &gr1, &ngr1)) != eslOK) cm_Fail(errbuf); 
  if((status = get_markup_from_msa(msa2, errbuf, &comment2, &ncomment2, &gf_tag2, &gf2, &ngf2, &gs_tag2, &gs2, &ngs2, &gr_tag2, &gr2, &ngr2)) != eslOK) cm_Fail(errbuf); 

  /* convert alignments to parsetrees, and free alignments */
  if((status = Alignment2Parsetrees(msa1, cm, mtr1, errbuf, &sq1, &tr1)) != eslOK) cm_Fail(errbuf);
  msa1_nseq = msa1->nseq;
  esl_msa_Destroy(msa1);

  if((status = Alignment2Parsetrees(msa2, cm, mtr2, errbuf, &sq2, &tr2)) != eslOK) cm_Fail(errbuf);
  msa2_nseq = msa2->nseq;
  esl_msa_Destroy(msa2);

  /* merge parsetrees */
  nmerged_tr = msa1_nseq;
  if((status = append_parsetrees(tr2, msa2_nseq, &(tr1), &(nmerged_tr), errbuf)) != eslOK) cm_Fail(errbuf);
  /* merge sequences */
  nmerged_sq = msa1_nseq;
  if((status = append_sequences(sq2, msa2_nseq, &(sq1), &(nmerged_sq), errbuf)) != eslOK) cm_Fail(errbuf);

  /* parsetrees to alignment */
  if((status = Parsetrees2Alignment(cm, abc_out, sq1, NULL, tr1, msa1_nseq + msa2_nseq, (! esl_opt_GetBoolean(go, "--resonly")), esl_opt_GetBoolean(go, "--matchonly"), &merged_msa)) != eslOK)
    cm_Fail("Error creating alignment from merged parsetrees.");

  /* add comments, GF, GS and GR (not GC) markup from initial alignments to merged alignment */
  if((status = add_msa_markup(merged_msa, errbuf, msa1_nseq, 0,         comment1, ncomment1, gf_tag1, gf1, ngf1, gs_tag1, gs1, ngs1, gr_tag1, gr1, ngr1)) != eslOK) cm_Fail(errbuf);
  if((status = add_msa_markup(merged_msa, errbuf, msa2_nseq, msa1_nseq, comment2, ncomment2, gf_tag2, gf2, ngf2, gs_tag2, gs2, ngs2, gr_tag2, gr2, ngr2)) != eslOK) cm_Fail(errbuf);

  /* output alignment */
  status = esl_msa_Write(ofp, merged_msa, (esl_opt_GetBoolean(go, "-1") ? eslMSAFILE_PFAM : eslMSAFILE_STOCKHOLM));
  if      (status == eslEMEM) cm_Fail("Memory error when outputting alignment\n");
  else if (status != eslOK)   cm_Fail("Writing alignment file failed with error %d\n", status);

  /* clean up */
  CMFileClose(cmfp);
  esl_msafile_Close(msafp1);
  esl_msafile_Close(msafp2);
  FreeParsetree(mtr1);
  FreeParsetree(mtr2);
  for(i = 0; i < nmerged_tr; i++) { 
    FreeParsetree(tr1[i]);
    esl_sq_Destroy(sq1[i]);
  }    
  esl_msa_Destroy(merged_msa);
  free(tr1);
  free(tr2);
  free(sq1);
  free(sq2);
  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_alphabet_Destroy(abc_out);
  return;
}

/* append_parsetrees()
 *
 * Given two arrays of parsetrees, swap pointers to append
 * the parsetrees in one array (<tr_to_append>) to the other, 
 * (<(*ret_tr)>) increasing the size of the first array. 
 */
static int
append_parsetrees(Parsetree_t **tr_to_append, int ntr_to_append, Parsetree_t ***ret_tr, int *ret_ntr, char *errbuf)
{
  int status;
  int ntr, ntr_existing, i, ip;
  void *tmp;
  ntr_existing = *ret_ntr;
  ntr = ntr_existing + ntr_to_append;
  
  
  ESL_RALLOC((*ret_tr), tmp, sizeof(Parsetree_t *) * ntr);
  for(i = ntr_existing; i < ntr; i++) { 
    ip = i - ntr_existing;
    (*ret_tr)[i] = tr_to_append[ip];
  }
  *ret_ntr = ntr;
  return eslOK;

 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "Memory reallocation error in append_parsetrees()");
  return status; /* NEVERREACHED */
}

/* append_sequences()
 *
 * Given two arrays of sequences, swap pointers to append
 * the sequences in one array (<sq_to_append>) to the other, 
 * (<(*ret_sq)>) increasing the size of the first array. 
 */
static int
append_sequences(ESL_SQ **sq_to_append, int nsq_to_append, ESL_SQ ***ret_sq, int *ret_nsq, char *errbuf)
{
  int status;
  int nsq, nsq_existing, i, ip;
  void *tmp;
  nsq_existing = *ret_nsq;
  nsq = nsq_existing + nsq_to_append;

  ESL_RALLOC((*ret_sq), tmp, sizeof(ESL_SQ *) * nsq);
  for(i = nsq_existing; i < nsq; i++) { 
    ip = i - nsq_existing;
    (*ret_sq)[i] = sq_to_append[ip];
  }
  *ret_nsq = nsq;
  return eslOK;

 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "Memory reallocation error in append_sequences()");
  return status; /* NEVERREACHED */
}

/* get_markup_from_msa
 *                   
 * Given an MSA return the markup that is unparsed by the esl_msa module 
 * (comments, GF, GS, GR, markup), with the exception of GC markup (although
 * parsed GC markup, RF and SS_cons IS included). The reason we don't include
 * GC markup is that there is no guarantee all columns of original alignments
 * will appear in the merged alignment because all gap columns will be removed
 * (this is due to the Alignment2Parsetrees(), Parsetrees2Alignment() strategy
 * we use for --merge). Also, it would be a pain to map each column of each
 * input msa to the merged msa (though doable, it requires looking at each
 * residue of each sequence). And finally, --merge is meant to merge two alignments
 * output by cmalign, which should not include any non-parsed #=GC markup,
 * so this isn't a huge deal.
 */
static int
get_markup_from_msa(ESL_MSA *msa, char *errbuf, char ***ret_comment, int *ret_ncomment, char ***ret_gf_tag, char ***ret_gf, int *ret_ngf, 
		    char ***ret_gs_tag, char ****ret_gs, int *ret_ngs, char ***ret_gr_tag, char ****ret_gr, int *ret_ngr)
{
  int status;

  int ncomment;
  char **comment;

  int   ngf;
  char **gf_tag;
  char **gf;

  int    ngs;
  char  **gs_tag;
  char ***gs;

  int    ngr;
  char  **gr_tag;
  char ***gr;

  int i,j;

  /* comment */
  if(msa->ncomment == 0) { 
    ncomment = 0;
    comment = NULL;
  }
  else { 
    ncomment = msa->ncomment;
    ESL_ALLOC(comment, sizeof(char *) * ncomment);
    for(i = 0; i < ncomment; i++) esl_strdup(msa->comment[i], -1, &(comment[i]));
  }

  /* gf */
  if(msa->ngf == 0) { 
    ngf = 0;
    gf_tag = NULL;
    gf     = NULL;
  }
  else { 
    ngf = msa->ngf;
    ESL_ALLOC(gf_tag, sizeof(char *) * ngf);
    ESL_ALLOC(gf,     sizeof(char *) * ngf);
    for(i = 0; i < ngf; i++) { 
      esl_strdup(msa->gf_tag[i], -1, &(gf_tag[i]));
      esl_strdup(msa->gf[i],     -1, &(gf[i]));
    }
  }

  /* gs */
  if(msa->ngs == 0) { 
    ngs = 0;
    gs_tag = NULL;
    gs     = NULL;
  }
  else { 
    ngs = msa->ngs;
    ESL_ALLOC(gs_tag, sizeof(char *)  * ngs);
    ESL_ALLOC(gs,     sizeof(char **) * ngs);
    for(i = 0; i < ngs; i++) { 
      esl_strdup(msa->gs_tag[i], -1, &(gs_tag[i]));
      ESL_ALLOC(gs[i], sizeof(char *) * msa->nseq);
      for(j = 0; j < msa->nseq; j++) { 
	esl_strdup(msa->gs[i][j], -1, &(gs[i][j]));
      }
    }
  }

  /* gr */
  if(msa->ngr == 0) { 
    ngr = 0;
    gr_tag = NULL;
    gr     = NULL;
  }
  else { 
    ngr = msa->ngr;
    ESL_ALLOC(gr_tag, sizeof(char *)  * ngr);
    ESL_ALLOC(gr,     sizeof(char **) * ngr);
    for(i = 0; i < ngr; i++) { 
      esl_strdup(msa->gr_tag[i], -1, &(gr_tag[i]));
      ESL_ALLOC(gr[i], sizeof(char *) * msa->nseq);
      for(j = 0; j < msa->nseq; j++) { 
	esl_strdup(msa->gr[i][j], -1, &(gr[i][j]));
	esl_strdealign(gr[i][j], gr[i][j], "-_.~", NULL);
      }
    }
  }

  *ret_ncomment = ncomment;
  *ret_comment  = comment;

  *ret_ngf      = ngf;
  *ret_gf_tag   = gf_tag;
  *ret_gf       = gf;

  *ret_ngs      = ngs;
  *ret_gs_tag   = gs_tag;
  *ret_gs       = gs;

  *ret_ngr      = ngr;
  *ret_gr_tag   = gr_tag;
  *ret_gr       = gr;

  return eslOK;
  
 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "Memory allocation error.");
  return status;
}

/* add_msa_markup
 *                   
 * Add markup to a MSA. Note: we free the
 * markup as soon as we've added it to the msa, 
 * because adding it to the msa means duplicating it.
 */
static int add_msa_markup(ESL_MSA *merged_msa, char *errbuf, int nseq, int seq_offset,
			    char **comment, int ncomment,
			    char **gf_tag, char **gf,  int ngf,
			    char **gs_tag, char ***gs, int ngs,
			    char **gr_tag, char ***gr, int ngr)
{ 
  int status;
  int i, m, ip;
  char *tmp_s;

  /* comments */
  for(m = 0; m < ncomment; m++) { 
    if(comment[m] != NULL) { 
      esl_msa_AddComment(merged_msa, comment[m]);
      free(comment[m]);
    }
  }
  if(comment != NULL) free(comment);

  /* GF */
  for(m = 0; m < ngf; m++) { 
    if(gf[m] != NULL) { 
      esl_msa_AddGF(merged_msa, gf_tag[m], gf[m]);
      free(gf[m]);
    }
    if(gf_tag[m] != NULL) free(gf_tag[m]);
  }
  if(gf_tag != NULL) free(gf_tag);
  if(gf != NULL)     free(gf);

  /* GS */
  for(m = 0; m < ngs; m++) { 
    for(i = 0; i < nseq; i++) { 
      ip = i + seq_offset;
      if(gs[m][i] != NULL) { 
	esl_msa_AddGS(merged_msa, gs_tag[m], ip, gs[m][i]);
	free(gs[m][i]);
      }
    }
    if(gs[m] != NULL) free(gs[m]);
    if(gs_tag[m] != NULL) free(gs_tag[m]);
  }    
  if(gs_tag != NULL) free(gs_tag);
  if(gs != NULL)     free(gs);
  
  /* GR */
  for(m = 0; m < ngr; m++) { 
    for(i = 0; i < nseq; i++) { 
      /* add gaps to full length of alignment */
      if(gr[m][i] != NULL) { 
	ip = i + seq_offset;
	if((status = gapize_string_to_fit_alignment(gr[m][i], merged_msa->aseq[ip], merged_msa->alen, '.', "-_.~", errbuf, &tmp_s)) != eslOK) return status;
	free(gr[m][i]);
	esl_msa_AppendGR(merged_msa, gr_tag[m], ip, tmp_s);
	free(tmp_s);
      }
    }
    if(gr[m] != NULL) free(gr[m]);
    if(gr_tag[m] != NULL) free(gr_tag[m]);
  }
  if(gr_tag != NULL) free(gr_tag);
  if(gr != NULL)     free(gr);

  return eslOK;
}


/* Function:  gapize_string_to_fit_alignment()
 * Synopsis:  Align a string so it fits into a msa (has aligned length of the msa)
 *            according to gaps in a reference aseq.
 * Incept:    EPN, Thu Sep 11 05:16:35 2008
 *
 * Purpose:   Align string <s> in place, by adding gap character <gapchar>
 *            whereever gaps in <aln_gapchars> appear in <aseq>. Return
 *            the new sequence in <ret_news>. Caller must free s.
 *            
 * Args:      s        - string to align
 *            aseq     - reference aligned sequence seq
 *            alen     - length of aseq (required)
 *            gapchar_to_add - the gap char to add to s
 *            aln_gapchars   - the possible gap characters in aseq
 *            errbuf   - for error messages
 *            ret_news - new string, gapized to alen
 *
 * Returns:   <eslOK> on success. <eslEINCOMPAT> if number of non gaps
 *            in aseq does not equal unaligned length of s.
 */
int
gapize_string_to_fit_alignment(char *s, const char *aseq, int alen, char gapchar_to_add, char *aln_gapchars, char *errbuf, char **ret_news)
{
  int status;
  int64_t uapos = 0;
  int64_t apos;
  char *news;
  int ualen;
  
  ualen = strlen(s);
  ESL_ALLOC(news, sizeof(char) * (alen+1));
  for (apos = 0; apos < alen; apos++) { 
    if(uapos > ualen) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "gapize_string_to_fit_alignment(), unaligned length of aligned template string (%" PRId64 ") not equal to unaligned string length: %d.", uapos, ualen);
    news[apos] = (strchr(aln_gapchars, aseq[apos]) == NULL) ? s[uapos++] : gapchar_to_add;
  }

  if(s[uapos] != '\0') ESL_FAIL(eslEINCOMPAT, errbuf, "gapize_string_to_fit_alignment(), unaligned length of aligned template string (%" PRId64 ") not equal to unaligned string length: %d.", uapos, ualen);
  news[alen] = '\0';
  *ret_news = news;

  return eslOK;

 ERROR: 
  ESL_FAIL(eslEMEM, errbuf, "Memory allocation error.");
  return status; /* NEVERREACHED */
}

#ifdef HAVE_MPI
/* determine_nseq_per_worker()
 * Given a CM, return the number of sequences we think we should send
 * to each worker (we don't know the number of sequences in the file).
 * The calculation is based on trying to get a worker to spend 
 * a specific amount of time MPI_WORKER_ALIGN_TARGET_SEC, a constant
 * from structs.h. 
 */
static int
determine_nseq_per_worker(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, int *ret_nseq_worker)
{
  /**ret_nseq_worker = 5;*/
  *ret_nseq_worker = 1;
  return eslOK;
}

/* add_worker_seqs_to_master
 * Add results (parstrees or CP9 traces, and possibly postcodes) from a
 * worker's seqs_to_aln object to a master seqs_to_aln object.
 */
static int
add_worker_seqs_to_master(seqs_to_aln_t *master_seqs, seqs_to_aln_t *worker_seqs, int offset)
{
  int x;

  if(worker_seqs->sq != NULL) cm_Fail("add_worker_seqs_to_master(), worker_seqs->sq non-NULL.");
  if(master_seqs->nseq < (offset + worker_seqs->nseq)) cm_Fail("add_worker_seqs_to_master(), master->nseq: %d, offset %d, worker->nseq: %d\n", master_seqs->nseq, offset, worker_seqs->nseq);

  if(worker_seqs->tr != NULL) {
    if(master_seqs->tr == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned parsetrees, master->tr is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(master_seqs->tr[x] == NULL); 
      master_seqs->tr[x] = worker_seqs->tr[(x-offset)];
    }
  }

  if(worker_seqs->cp9_tr != NULL) {
    if(master_seqs->cp9_tr == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned cp9 traces, master->cp9_tr is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(master_seqs->cp9_tr[x] == NULL); 
      master_seqs->cp9_tr[x] = worker_seqs->cp9_tr[(x-offset)];
    }
  }

  if(worker_seqs->postcode1 != NULL) {
    if(master_seqs->postcode1 == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned postcodes, master->postcode1 is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(master_seqs->postcode1[x] == NULL); 
      master_seqs->postcode1[x] = worker_seqs->postcode1[(x-offset)];
    }
  }

  if(worker_seqs->postcode2 != NULL) {
    if(master_seqs->postcode2 == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned postcodes, master->postcode2 is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(master_seqs->postcode2[x] == NULL); 
      master_seqs->postcode2[x] = worker_seqs->postcode2[(x-offset)];
    }
  }

  if(worker_seqs->sc != NULL) {
    if(master_seqs->sc == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned scores, master->sc is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(!(NOT_IMPOSSIBLE(master_seqs->sc[x])));
      master_seqs->sc[x] = worker_seqs->sc[(x-offset)];
    }
  }

  if(worker_seqs->pp != NULL) {
    if(master_seqs->pp == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned post probs, master->pp is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(!(NOT_IMPOSSIBLE(master_seqs->pp[x])));
      master_seqs->pp[x] = worker_seqs->pp[(x-offset)];
    }
  }

  if(worker_seqs->struct_sc != NULL) {
    if(master_seqs->struct_sc == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned post probs, master->struct_sc is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(!(NOT_IMPOSSIBLE(master_seqs->struct_sc[x])));
      master_seqs->struct_sc[x] = worker_seqs->struct_sc[(x-offset)];
    }
  }

  return eslOK;
}
#endif /* of #ifdef HAVE_MPI */

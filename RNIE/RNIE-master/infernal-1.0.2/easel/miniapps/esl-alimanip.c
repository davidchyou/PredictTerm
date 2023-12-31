/* Manipulate a multiple sequence alignment in some useful ways.
 *
 * Derived from easel's esl-alistat which was from squid's alistat (1995)
 * EPN, Fri Aug 10 08:52:30 2007
 * SVN $Id: esl-alimanip.c 416 2009-10-22 11:42:55Z nawrockie $
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "easel.h"
#include "esl_distance.h"
#include "esl_fileparser.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_msa.h"
#include "esl_distance.h"
#include "esl_dmatrix.h"
#include "esl_vectorops.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_wuss.h"

static char banner[] = "manipulate a multiple sequence alignment file";
static char usage[]  = "[options] <msafile>\n\
The <msafile> must be in Stockholm format.";

#define OTHERMSAOPTS  "--merge,--morph"           /* Exclusive choice for options involving other MSAs */

static int  keep_or_remove_rf_gaps(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, int keep_flag, int remove_flag);
static int  write_rf_gapthresh(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa);
static int  write_rf_given_alen(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, char *amask, int amask_len);
static int  write_rf_given_rflen(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, char *rfmask, int rfmask_len);
static int  write_rf_given_useme(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, int *useme);
static int  morph_msa(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa1, ESL_MSA *msa2, ESL_MSA **newmsa1);
static int  merge_msa(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa1, ESL_MSA *msa2, ESL_MSA **ret_merged_msa);
static int  add_gap_columns_to_msa(char *errbuf, ESL_MSA *msa, int *toadd, ESL_MSA **ret_msa, int do_treat_as_rf_gap);
static int  cp_and_add_gaps_to_aseq(char *new_aseq, char *orig_aseq, int alen, int *toadd, int nnew, char gapchar);
static int  is_flush_left(int *ngaps, int astart, int aend);
static int  is_flush_right(int *ngaps, int astart, int aend);
static int  pick_gappiest_columns(int *ngaps, int astart, int aend, int nkeep, int **ret_cols_to_keep);
static int  get_gaps_per_column(ESL_MSA *msa, int **ret_ngaps);
static int  map_cpos_to_apos(ESL_MSA *msa, int **ret_c2a_map, int *ret_clen);
static int  individualize_consensus(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa);
static int  read_sqfile(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc, int nseq, ESL_SQ ***ret_sq);
static int  trim_msa(ESL_MSA *msa, ESL_SQ **sq, char *errbuf);
static int  dump_insert_info(FILE *fp, ESL_MSA *msa, char *errbuf);
static int  dump_residue_info(FILE *fp, ESL_MSA *msa, char *errbuf);
static int  dump_delete_info(FILE *fp, ESL_MSA *msa, char *errbuf);
static int  plot_inserts(FILE *fp, ESL_MSA *msa, int do_log, char *errbuf);
static int  dump_infocontent(FILE *fp, ESL_MSA *msa, char *errbuf);
static int  plot_gaps(FILE *fp, ESL_MSA *msa, char *errbuf);
static int  get_tree_order(ESL_TREE *T, char *errbuf, int **ret_order);
static int  reorder_msa(ESL_MSA *msa, int *order, char *errbuf);
static int  dmx_Visualize(FILE *fp, ESL_DMATRIX *D, double min, double max);
static int  read_mask_file(char *filename, char *errbuf, char **ret_mask, int *ret_mask_len);
static int  handle_post_opts(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa);
static int  output_rf_as_mask(FILE *fp, char *errbuf, ESL_MSA *msa);
static int  expand_msa2mask(char *errbuf, ESL_MSA *msa1, char *xmask, ESL_MSA **newmsa1);
static int  compare_ints(const void *el1, const void *el2);
static int  msa_median_length(ESL_MSA *msa);
static int  msa_remove_seqs_below_minlen(ESL_MSA *msa, float minlen, ESL_MSA **ret_new_msa);
static int  msa_remove_truncated_seqs(ESL_MSA *msa, char *errbuf, int ntrunc, ESL_MSA **ret_new_msa);
static int  number_columns(ESL_MSA *msa, int do_all, char *errbuf);
static char digit_to_char(int digit);
static int  int_ndigits(int i);
static char get_char_digit_x_from_int(int i, int place);
static int  keep_contiguous_column_block(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa);
static int  read_seq_name_file(char *filename, char *errbuf, char ***ret_seqlist, int *ret_seqlist_n);
static int  msa_keep_or_remove_seqs(ESL_MSA *msa, char *errbuf, char **seqlist, int seqlist_n, int do_keep, ESL_MSA **ret_new_msa);

static ESL_OPTIONS options[] = {
  /* name          type        default  env   range      togs reqs  incomp                      help                                                       docgroup */
  { "-h",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "help; show brief info on version and usage",                     1 },
  { "-o",          eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, NULL,                       "output the alignment to file <f>, not stdout",                   1 },
  { "-1",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "output alignment in Pfam (non-interleaved, 1 line/seq) format",  1 },
  { "--list",      eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, NULL,                       "output list of sequence names in alignment to file <f>",         1 },
  { "--devhelp",   eslARG_NONE,  NULL,  NULL, NULL,      NULL,NULL, NULL,                       "show list of undocumented developer options",                    1 },
  { "-g",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "add/rewrite #=GC RF markup based on gap frequency in each col",  2 },
  { "--gapthresh", eslARG_REAL,  "0.5", NULL, "0<=x<=1", NULL,"-g", NULL,                       "with -g, fraction of gaps allowed in non-gap RF columns [0.5]",  2 },
  { "--mask-all",  eslARG_INFILE,FALSE,NULL, NULL,      NULL,NULL, NULL,                        "set #=GC RF as x=1, gap=0 from 1/0s in 1-line <f> (len=alen)",   2 },
  { "--mask-rf",   eslARG_INFILE, FALSE,NULL, NULL,      NULL,NULL, NULL,                       "set #=GC RF as x=1, gap=0 from 1/0s in 1-line <f> (len=rf len)", 2 },
  { "--pfract",    eslARG_REAL,  NULL,  NULL, "0<=x<=1", NULL,NULL, NULL,                       "set #=GC RF as cols w/<x> fraction of seqs w/POST >= --pthresh", 2 },
  { "--pthresh",   eslARG_REAL,  "0.9", NULL, "0<=x<=1", NULL,"--pfract", NULL,                 "set #=GR POST threshold for --pfract as <x> [default=0.9]",      2 },
  { "--p-rf",      eslARG_NONE,  NULL,  NULL, NULL,      NULL,"--pfract", NULL,                 "with --pfract options, ignore gap #=GC RF columns",              2 },
  { "-k",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "keep  only non-gap RF columns, as RF is defined in input aln", 3 },
  { "-r",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "remove all non-gap RF columns, as RF is defined in input aln", 3 },
  { "--kmask",     eslARG_OUTFILE,FALSE,NULL, NULL,      NULL,"-k", NULL,                       "w/-k, output RF as mask to <f> before removing gap columns", 3},
  { "--start-all", eslARG_INT,   NULL,  NULL, NULL,      NULL,"--end-all",  "--start-rf",       "keep columns starting at column <n>", 3 },
  { "--end-all",   eslARG_INT,   NULL,  NULL, NULL,      NULL,"--start-all","--start-rf",       "keep columns ending   at column <n>", 3 },
  { "--start-rf",  eslARG_INT,   NULL,  NULL, NULL,      NULL,"--end-rf",   "--start-all",      "keep columns starting at non-gap RF column <n>", 3 },
  { "--end-rf",    eslARG_INT,   NULL,  NULL, NULL,      NULL,"--start-rf", "--start-all",      "keep columns ending   at non-gap RF column <n>", 3 },
  { "--tree",      eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL,OTHERMSAOPTS,                "reorder MSA to tree order following SLC, save Newick tree to <f>", 4 },
  { "--lfract",    eslARG_REAL,  NULL,  NULL, "0<=x<=1", NULL,NULL, NULL,                       "remove sequences w/length < <x> fraction of median length",      4 },
  { "--lmin",      eslARG_INT,   NULL,  NULL, "n>0",     NULL,NULL, NULL,                       "remove sequences w/length < <n> residues",                       4 },
  { "--detrunc",   eslARG_INT,   NULL,  NULL, "n>0",     NULL,NULL, NULL,                       "remove seqs w/gaps in >= <n> 5' or 3'-most non-gap #=GC RF cols",4 },
  { "--seq-r",     eslARG_INFILE,NULL,  NULL, NULL,      NULL,NULL, "--seq-k",                  "remove sequences with names listed in file <f>",                 4 },
  { "--seq-k",     eslARG_INFILE,NULL,  NULL, NULL,      NULL,NULL, "--seq-r",                  "remove all sequences *except* those listed in file <f>",         4 },
  { "--trim",      eslARG_INFILE, NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "trim aligned seqs in <msafile> to subseqs in <f>",               4 },
  { "--iinfo",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,                "print info on # of insertions b/t all non-gap RF cols to <f>",   5 },
  { "--icinfo",    eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "print info on information content of each non-gap RF column",    5 },
  { "--rinfo",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "print info on # of residues in each col of alignment to <f>",    5 },
  { "--dinfo",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "print info on # of deletes in non-gap RF cols of aln to <f>",    5 },
  { "--pinfo",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "print info on posterior probabilities in <msafile> to <f>",      5 },
  { "--sindi",     eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, "-g,-k,-r,--morph",         "annotate individual secondary structures by imposing consensus", 7 },
  { "--num-all",   eslARG_NONE,   NULL, NULL, NULL,      NULL,NULL, NULL,                       "add annotation numbering all columns",                          11 },
  { "--num-rf",    eslARG_NONE,   NULL, NULL, NULL,      NULL,NULL, NULL,                       "add annotation numbering the non-gap RF columns",               11 },
  { "--omask",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, NULL,                       "output RF annotation as 1/0 mask to file <f>",                   9 },
  { "--amino",     eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL,"--dna,--rna",               "<msafile> contains protein alignments",                         10 },
  { "--dna",       eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL,"--amino,--rna",             "<msafile> contains DNA alignments",                             10 },
  { "--rna",       eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL,"--amino,--dna",             "<msafile> contains RNA alignments",                             10 },

  /* All options below are developer options, only shown if --devhelp invoked */
  { "--iplot",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL,OTHERMSAOPTS,                "plot heatmap of # of insertions b/t all non-gap RF cols to <f>", 101 },
  { "--ilog",      eslARG_NONE,  FALSE, NULL, NULL,      NULL,"--iplot", NULL,                  "w/--iplot, use log scale for heatmap of insert counts",          101 },
  { "--gplot",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "plot checkerboard grid of # of gaps in non-gap RF cols to <f>",  101 },
  { "--morph",     eslARG_INFILE, NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "morph msa in <msafile> to msa in <f>'s gap structure",          101 },
  { "--merge",     eslARG_INFILE,FALSE, NULL, NULL,      NULL,NULL, "--morph,-g,-k,-r",         "merge msa in <msafile> with msa in <f>",                         101 },

  { "--xmask",     eslARG_INFILE, NULL, NULL, NULL,      NULL,NULL, NULL,                       "for each 0 column in <f>, add a 100% gap column to <msafile>",   102 },
  { "--verbose",   eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "be verbose (usually with --morph or --merge)",            102 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = NULL;	/* application configuration       */
  ESL_ALPHABET *abc     = NULL;	/* biological alphabet             */
  char         *alifile = NULL;	/* alignment file name             */
  int           fmt;		/* format code for alifile         */
  ESL_MSAFILE  *afp     = NULL;	/* open alignment file             */
  ESL_MSA      *msa     = NULL;	/* one multiple sequence alignment */
  int           status;		/* easel return code               */
  int           nali;		/* number of alignments read       */
  FILE         *ofp;		/* output file (default is stdout) */
  char          errbuf[eslERRBUFSIZE];
  int           write_ali = FALSE; /* set to TRUE if we should print a new MSA */
  /* --merge, --morph related vars */
  ESL_MSAFILE  *otherafp = NULL;	/* other input alignment file (with --morph) */
  ESL_MSA      *othermsa = NULL;	/* other input alignment      (with --morph) */
  /* --trim related vars */
  ESL_SQFILE   *trimfp = NULL;  /* sequence file with subsequences for --trim */
  /* --iinfo, --iplot, --gplot --rinfo, --dinfo related vars */
  FILE *treefp  = NULL;  /* output file for --tree */
  FILE *iinfofp = NULL;  /* output file for --iinfo */
  FILE *iplotfp = NULL;  /* output file for --iplot */
  FILE *gplotfp = NULL;  /* output file for --gplot */
  FILE *rinfofp = NULL;  /* output file for --rinfo */
  FILE *dinfofp = NULL;  /* output file for --dinfo */
  FILE *icinfofp = NULL; /* output file for --icinfo */
  FILE *listfp = NULL;   /* output file for --list */
  /* --mask-all */
  char *amask = NULL;
  int   amask_len = -1;
  /* --mask-all */
  char *rfmask = NULL;
  int   rfmask_len = -1;
  /* --xmask */
  char *xmask = NULL;
  int   xmask_len = -1;
  /* --omask */
  FILE *omaskfp;
  /* --kmask */
  FILE *kmaskfp;

  /***********************************************
   * Parse command line
   ***********************************************/

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  if (esl_opt_GetBoolean(go, "--devhelp") )
    {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      puts("\noptions for adding/rewriting #=GC RF annotation:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\noptions for removing columns:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
      puts("\noptions for numbering columns:");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 80); 
      puts("\noptions for reordering/removing/trimming sequences:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
      puts("\noptions for displaying info on inserts/gaps/posterior probabilities:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 
      puts("\noptions for manipulating secondary structure annotation:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 
      puts("\noptions for outputting a lanemask file:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
      puts("\noptions for specifying input alphabet:");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80);
      puts("\nundocumented, experimental developer options:");
      esl_opt_DisplayHelp(stdout, go, 101, 2, 80);
      puts("\noptions for comparison/modification based on another MSA file:");
      esl_opt_DisplayHelp(stdout, go, 102, 2, 80); 
      exit(0);
    }
  if (esl_opt_GetBoolean(go, "-h") )
    {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      puts("\noptions for adding/rewriting #=GC RF annotation:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\noptions for removing columns:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
      puts("\noptions for numbering columns:");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 80); 
      puts("\noptions for reordering/removing/trimming sequences:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
      puts("\noptions for displaying info on inserts/gaps/posterior probabilities:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 
      puts("\noptions for manipulating secondary structure annotation:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 
      puts("\noptions for outputting a lanemask file:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
      puts("\noptions for specifying input alphabet:");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80);
      exit(0);
    }

  if (esl_opt_ArgNumber(go) != 1) 
    {
      printf("Incorrect number of command line arguments.\n");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  alifile = esl_opt_GetArg(go, 1);

  fmt             = eslMSAFILE_STOCKHOLM;

  /***********************************************
   * Open the MSA file; determine alphabet; set for digital input
   ***********************************************/

  status = esl_msafile_Open(alifile, fmt, NULL, &afp);
  if (status == eslENOTFOUND) 
    esl_fatal("Alignment file %s doesn't exist or is not readable\n", alifile);
  else if (status == eslEFORMAT) 
    esl_fatal("Couldn't determine format of alignment %s\n", alifile);
  else if (status != eslOK) 
    esl_fatal("Alignment file open failed with error %d\n", status);

  /* open output file */
  if (esl_opt_GetString(go, "-o") != NULL) {
    if ((ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
    } else ofp = stdout;

  if      (esl_opt_GetBoolean(go, "--amino"))   abc = esl_alphabet_Create(eslAMINO);
  else if (esl_opt_GetBoolean(go, "--dna"))     abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--rna"))     abc = esl_alphabet_Create(eslRNA);
  else {
    int type;
    status = esl_msafile_GuessAlphabet(afp, &type);
    if (status == eslEAMBIGUOUS)    esl_fatal("Failed to guess the bio alphabet used in %s.\nUse --dna, --rna, or --amino option to specify it.", alifile);
    else if (status == eslEFORMAT)  esl_fatal("Alignment file parse failed: %s\n", afp->errbuf);
    else if (status == eslENODATA)  esl_fatal("Alignment file %s is empty\n", alifile);
    else if (status != eslOK)       esl_fatal("Failed to read alignment file %s\n", alifile);
    abc = esl_alphabet_Create(type);
  }
  esl_msafile_SetDigital(afp, abc);
  if((esl_opt_GetBoolean(go, "--sindi")) && (abc->type != eslRNA && abc->type != eslDNA))
    esl_fatal("-i option pertains to base pairs and only makes sense with DNA or RNA alphabets.");

  /* optionally, open --morph or --merge msa file for reading, --merge and --morph are incompatible
   * with each other, so we'll never try to do open othermsafile more than once.
   */
  if(esl_opt_GetString(go, "--morph") != NULL)
    {
      status = esl_msafile_OpenDigital(abc, esl_opt_GetString(go, "--morph"), eslMSAFILE_STOCKHOLM, NULL, &otherafp);
      if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "--morph alignment file %s doesn't exist or is not readable\n", 
					      esl_opt_GetString(go, "--morph"));
      else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of --morph alignment %s\n", 
					      esl_opt_GetString(go, "--morph"));
      else if (status != eslOK)      ESL_FAIL(status, errbuf, "Alignment file open failed with error %d\n", status);
    }
  if(esl_opt_GetString(go, "--merge") != NULL)
    {
      status = esl_msafile_OpenDigital(abc, esl_opt_GetString(go, "--merge"), eslMSAFILE_STOCKHOLM, NULL, &otherafp);
      if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "--merge alignment file %s doesn't exist or is not readable\n", 
					      esl_opt_GetString(go, "--merge"));
      else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of --merge alignment %s\n", 
					      esl_opt_GetString(go, "--merge"));
      else if (status != eslOK)      ESL_FAIL(status, errbuf, "Alignment file open failed with error %d\n", status);
    }

  /* read --mask-all file, if nec */
  if(esl_opt_GetString(go, "--mask-all") != NULL) {
    if((status = read_mask_file(esl_opt_GetString(go, "--mask-all"), errbuf, &amask, &amask_len)) != eslOK)
      esl_fatal(errbuf);
  }
  /* read --mask-rf file, if nec */
  if(esl_opt_GetString(go, "--mask-rf") != NULL) {
    if((status = read_mask_file(esl_opt_GetString(go, "--mask-rf"), errbuf, &rfmask, &rfmask_len)) != eslOK)
      esl_fatal(errbuf);
  }
  /* read --xmask file, if nec */
  if(esl_opt_GetString(go, "--xmask") != NULL) {
    if((status = read_mask_file(esl_opt_GetString(go, "--xmask"), errbuf, &xmask, &xmask_len)) != eslOK)
      esl_fatal(errbuf);
  }
  /***********************************************
   * Read MSAs one at a time.
   ***********************************************/

  nali = 0;
  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;

      /* handle the --lfract option if enabled, all subsequent manipulations will omit any short seqs removed here */
      if(! esl_opt_IsDefault(go, "--lfract")) {
	int median   = msa_median_length(msa);
	float minlen = esl_opt_GetReal(go, "--lfract") * (float) median;
	ESL_MSA *new_msa;
	msa_remove_seqs_below_minlen(msa, minlen, &new_msa);
	/* new_msa is msa without seqs below minlen, swap ptrs */
	esl_msa_Destroy(msa);
	msa = new_msa;
	write_ali = TRUE;
      }


      /* handle the --lmin option if enabled, all subsequent manipulations will omit any short seqs removed here */
      if(! esl_opt_IsDefault(go, "--lmin")) {
	float minlen = esl_opt_GetInteger(go, "--lmin");
	ESL_MSA *new_msa;
	msa_remove_seqs_below_minlen(msa, minlen, &new_msa);
	/* new_msa is msa without seqs below minlen, swap ptrs */
	esl_msa_Destroy(msa);
	msa = new_msa;
	write_ali = TRUE;
      }

      /* handle the --detrunc option if enabled, all subsequent manipulations will omit any seqs removed here */
      if(! esl_opt_IsDefault(go, "--detrunc")) {
	ESL_MSA *new_msa;
	if((status =msa_remove_truncated_seqs(msa, errbuf, esl_opt_GetInteger(go, "--detrunc"), &new_msa)) != eslOK) esl_fatal(errbuf);
	/* new_msa is msa without seqs below minlen, swap ptrs */
	esl_msa_Destroy(msa);
	msa = new_msa;
	write_ali = TRUE;
      }

      /* handle the --seq-k and --seq-r options if enabled, all subsequent manipulations will omit any seqs removed here */
      if((! esl_opt_IsDefault(go, "--seq-k")) || (! esl_opt_IsDefault(go, "--seq-r"))) {
	ESL_MSA *new_msa;
	char   **seqlist;
	int      seqlist_n, n;
	if(! esl_opt_IsDefault(go, "--seq-k")) { 
	  if((status = read_seq_name_file(esl_opt_GetString(go, "--seq-k"), errbuf, &seqlist, &seqlist_n)) != eslOK) esl_fatal(errbuf);	  
	  if((status = msa_keep_or_remove_seqs(msa, errbuf, seqlist, seqlist_n, TRUE, &new_msa)) != eslOK)        esl_fatal(errbuf);	  
	  /* new_msa is msa but only with seqs listed in --seq-k <f> file */
	}
	else { /* --seq-r enabled */
	  if((status = read_seq_name_file(esl_opt_GetString(go, "--seq-r"), errbuf, &seqlist, &seqlist_n)) != eslOK) esl_fatal(errbuf);	  
	  if((status = msa_keep_or_remove_seqs(msa, errbuf, seqlist, seqlist_n, FALSE, &new_msa)) != eslOK)        esl_fatal(errbuf);	  
	  /* new_msa is msa but without seqs listed in --seq-r <f> file */
	}
	esl_msa_Destroy(msa);
	msa = new_msa;
	for(n = 0; n < seqlist_n; n++) free(seqlist[n]); 
	free(seqlist);
	write_ali = TRUE;
      }

      /* read other msa if --morph or --merge (which are incompatible with each other) is enabled */
      if((esl_opt_GetString(go, "--morph") != NULL) || (esl_opt_GetString(go, "--merge") != NULL))
	{
	  if ((status = esl_msa_Read(otherafp, &othermsa)) != eslOK) {
	    if(status == eslEFORMAT) 
	      esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", 
			otherafp->linenumber, otherafp->fname, otherafp->errbuf, otherafp->buf);	
	    else if (status == eslEOF)
	      esl_fatal("No alignments read in %s.", esl_opt_GetString(go, "--morph"));
	  }
	}

      /* if nec, handle --trim option */
      if(esl_opt_GetString(go, "--trim") != NULL) { 
	/* open seq file for --trim */
	status = esl_sqfile_Open(esl_opt_GetString(go, "--trim"), eslSQFILE_UNKNOWN, NULL, &(trimfp));
	if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "File %s doesn't exist or is not readable\n", esl_opt_GetString(go, "--trim"));
	else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of sequence file %s\n", esl_opt_GetString(go, "--trim"));
	else if (status == eslEINVAL)  ESL_FAIL(status, errbuf, "Can’t autodetect stdin or .gz."); 
	else if (status != eslOK)      ESL_FAIL(status, errbuf, "Sequence file open failed with error %d\n", status);
	/* read the sequences */
	ESL_SQ **sq;
	read_sqfile(trimfp, msa->abc, msa->nseq, &sq); /* dies on failure */
	/* trim the msa */
	if((status = trim_msa(msa, sq, errbuf)) != eslOK) goto ERROR;
	write_ali = TRUE;
      }

      /* if nec, morph <msafile> into gap structure in <f> (from --morph <f>)*/
      if(esl_opt_GetString(go, "--morph") != NULL)
	{
	  ESL_MSA *newmsa;
	  if((status = morph_msa(go, errbuf, msa, othermsa, &newmsa)) != eslOK) goto ERROR;
	  write_ali = TRUE;
	  /*status = esl_msa_Write(stdout, othermsa, eslMSAFILE_STOCKHOLM);
	    status = esl_msa_Write(stdout, newmsa, eslMSAFILE_STOCKHOLM);*/
	  msa = newmsa;
	}

      /* if nec, merge <msafile> and <f> (from --merge <f>) */
      if(esl_opt_GetString(go, "--merge") != NULL)
	{
	  ESL_MSA *newmsa;
	  if((status = merge_msa(go, errbuf, msa, othermsa, &newmsa)) != eslOK) goto ERROR;
	  write_ali = TRUE;
	  /*status = esl_msa_Write(stdout, othermsa, eslMSAFILE_STOCKHOLM);
	    status = esl_msa_Write(stdout, newmsa, eslMSAFILE_STOCKHOLM);*/
	  msa = newmsa;
	}

      /* rewrite RF annotation, if nec */
      if(esl_opt_GetBoolean(go, "-g")) {
	if((status = write_rf_gapthresh(go, errbuf, msa)) != eslOK) goto ERROR;
	write_ali = TRUE;
      }
      if(amask != NULL) { /* --mask-all enabled */
	if((status = write_rf_given_alen(go, errbuf, msa, amask, amask_len)) != eslOK) goto ERROR;
	write_ali = TRUE;
      }
      if(rfmask != NULL) { /* --mask-rf enabled */
	if((status = write_rf_given_rflen(go, errbuf, msa, rfmask, rfmask_len)) != eslOK) goto ERROR;
	write_ali = TRUE;
      }

      /* handle posterior (--p*) options, if nec */
      if(! ((esl_opt_IsDefault(go, "--pfract")) || (! (esl_opt_IsDefault(go, "--pinfo"))))) { 
	if((status = handle_post_opts(go, errbuf, msa) != eslOK)) goto ERROR;
	if(! (esl_opt_IsDefault(go, "--pfract")))
	  write_ali = TRUE;
      }

      /* Remove columns based on --start-all --end-all, --start-rf --end-rf, if nec */
      if((! esl_opt_IsDefault(go, "--start-all")) || (! esl_opt_IsDefault(go, "--start-rf")))
	{
	  if((status = keep_contiguous_column_block(go, errbuf, msa) != eslOK)) goto ERROR;
	  write_ali = TRUE;
	}

      /* keep or remove columns based on RF annotation, if nec */
      if(esl_opt_GetBoolean(go, "-k") || esl_opt_GetBoolean(go, "-r"))
	{
	  if(esl_opt_GetString(go, "--kmask") != NULL)
	    {
	      if ((kmaskfp = fopen(esl_opt_GetString(go, "--kmask"), "w")) == NULL) 
		ESL_FAIL(eslFAIL, errbuf, "Failed to open --kmask output file %s\n", esl_opt_GetString(go, "--kmask"));
	      if((status = output_rf_as_mask(kmaskfp, errbuf, msa)) != eslOK) goto ERROR;
	    }
	  if((status = keep_or_remove_rf_gaps(go, errbuf, msa, 
					      esl_opt_GetBoolean(go, "-k"),
					      esl_opt_GetBoolean(go, "-r"))) != eslOK) goto ERROR;
	  write_ali = TRUE;
	}

      /* impose consensus structure to get individual secondary structures, if nec */
      if(esl_opt_GetBoolean(go, "--sindi"))
	{
	  if((status = individualize_consensus(go, errbuf, msa) != eslOK)) goto ERROR;
	  write_ali = TRUE;
	}

      /* handle the --tree option, if enabled */
      if(! esl_opt_IsDefault(go, "--tree"))
	{
	  if ((treefp = fopen(esl_opt_GetString(go, "--tree"), "w")) == NULL) 
	    ESL_FAIL(eslFAIL, errbuf, "Failed to open --tree output file %s\n", esl_opt_GetString(go, "--tree"));

	  ESL_TREE    *T = NULL;/* the tree, created by Single-Linkage Clustering */
	  ESL_DMATRIX *D = NULL;/* the distance matrix */
	  
	  /* Create distance matrix and infer tree by single linkage clustering */
	  esl_dst_XDiffMx(msa->abc, msa->ax, msa->nseq, &D);
	  esl_tree_SingleLinkage(D, &T);
	  esl_tree_SetTaxaParents(T);
	  esl_tree_SetTaxonlabels(T, msa->sqname);

	  esl_tree_WriteNewick(treefp, T); 
	  fclose(treefp);
	  /*printf("# Tree saved in Newick format to file %s.\n", esl_opt_GetString(go, "--tree")); */

	  esl_tree_Validate(T, NULL);

	  /* Get new order for seqs in the MSA based on the tree */
	  int *order;
	  if((status = get_tree_order(T, errbuf, &order)) != eslOK) goto ERROR;
	  /*for(i = 0; i < msa->nseq; i++) {
	    printf("new MSA idx: %3d | orig MSA idx: %3d\n", i, order[i]);
	    }*/
	  esl_tree_Destroy(T);
	  esl_dmatrix_Destroy(D);
	  if((status = reorder_msa(msa, order, errbuf)) != eslOK) goto ERROR;
	  write_ali = TRUE;
	  free(order);
	}	  

      /* --xmask option: expand the alignment to fit lanemask in xmask <f>, number of TOTAL msa
       * columns must equal number of 1s in <f>.
       */
      if(xmask != NULL) { 
	ESL_MSA *newmsa;
	if((status = expand_msa2mask(errbuf, msa, xmask, &newmsa)) != eslOK) goto ERROR;
	  write_ali = TRUE;
	  msa = newmsa;
      }

      /* handle the --iinfo option, if enabled, do this after all MSA has been manipulated due to other options */
      if(! esl_opt_IsDefault(go, "--iinfo")) {
	if ((iinfofp = fopen(esl_opt_GetString(go, "--iinfo"), "w")) == NULL) 
	  ESL_FAIL(eslFAIL, errbuf, "Failed to open --iinfo output file %s\n", esl_opt_GetString(go, "--iinfo"));
	if((status = dump_insert_info(iinfofp, msa, errbuf) != eslOK)) goto ERROR;
	/*printf("# Insert information saved to file %s.\n", esl_opt_GetString(go, "--iinfo")); */
	fclose(iinfofp);
      }

      /* handle the --iplot option, if enabled, do this after all MSA has been manipulated due to other options */
      if(! esl_opt_IsDefault(go, "--iplot")) {
	if ((iplotfp = fopen(esl_opt_GetString(go, "--iplot"), "w")) == NULL) 
	  ESL_FAIL(eslFAIL, errbuf, "Failed to open --iplot output file %s\n", esl_opt_GetString(go, "--iplot"));
	if((status = plot_inserts(iplotfp, msa, esl_opt_GetBoolean(go, "--ilog"), errbuf) != eslOK)) goto ERROR;
	fclose(iplotfp);
      }

      /* handle the --icinfo option, if enabled, do this after all MSA has been manipulated due to other options */
      if(! esl_opt_IsDefault(go, "--icinfo")) {
	if ((icinfofp = fopen(esl_opt_GetString(go, "--icinfo"), "w")) == NULL) 
	  ESL_FAIL(eslFAIL, errbuf, "Failed to open --icinfo output file %s\n", esl_opt_GetString(go, "--iplot"));
	if((status = dump_infocontent(icinfofp, msa, errbuf) != eslOK)) goto ERROR;
	fclose(icinfofp);
      }

      /* handle the --gplot option, if enabled, do this after all MSA has been manipulated due to other options */
      if(! esl_opt_IsDefault(go, "--gplot")) {
	if ((gplotfp = fopen(esl_opt_GetString(go, "--gplot"), "w")) == NULL) 
	  ESL_FAIL(eslFAIL, errbuf, "Failed to open --gplot output file %s\n", esl_opt_GetString(go, "--gplot"));
	if((status = plot_gaps(gplotfp, msa, errbuf) != eslOK)) goto ERROR;
	fclose(gplotfp);
      }

      /* handle the --rinfo option, if enabled, do this after all MSA has been manipulated due to other options */
      if(! esl_opt_IsDefault(go, "--rinfo")) {
	if ((rinfofp = fopen(esl_opt_GetString(go, "--rinfo"), "w")) == NULL) 
	  ESL_FAIL(eslFAIL, errbuf, "Failed to open --rinfo output file %s\n", esl_opt_GetString(go, "--rinfo"));
	if((status = dump_residue_info(rinfofp, msa, errbuf) != eslOK)) goto ERROR;
	fclose(rinfofp);
      }

      /* handle the --dinfo option, if enabled, do this after all MSA has been manipulated due to other options */
      if(! esl_opt_IsDefault(go, "--dinfo")) {
	if ((dinfofp = fopen(esl_opt_GetString(go, "--dinfo"), "w")) == NULL) 
	  ESL_FAIL(eslFAIL, errbuf, "Failed to open --dinfo output file %s\n", esl_opt_GetString(go, "--dinfo"));
	if((status = dump_delete_info(dinfofp, msa, errbuf) != eslOK)) goto ERROR;
	fclose(dinfofp);
      }

      /* handle the --num-rf and --num-all options, if enabled, do this after all MSA has been manipulated due to other options */
      if(! esl_opt_IsDefault(go, "--num-rf")) { 
	if((status = number_columns(msa, FALSE, errbuf) != eslOK)) goto ERROR;
	write_ali = TRUE;
      }
      if(! esl_opt_IsDefault(go, "--num-all")) { 
	if((status = number_columns(msa, TRUE, errbuf) != eslOK)) goto ERROR;
	write_ali = TRUE;
      }

      /* write out list of sequences, if nec */
      if(! esl_opt_IsDefault(go, "--list")) {
	if ((listfp = fopen(esl_opt_GetString(go, "--list"), "w")) == NULL) 
	  esl_fatal("Failed to open --list output file %s\n", esl_opt_GetString(go, "--list"));
	int i;
	for(i = 0; i < msa->nseq; i++) fprintf(listfp, "%s\n", msa->sqname[i]);
	fclose(listfp);
     } 

      /* write out alignment, if nec */
      if(write_ali || esl_opt_GetBoolean(go, "-1")) {
	status = esl_msa_Write(ofp, msa, (esl_opt_GetBoolean(go, "-1") ? eslMSAFILE_PFAM : eslMSAFILE_STOCKHOLM));
	if      (status == eslEMEM) ESL_FAIL(status, errbuf, "Memory error when outputting alignment\n");
	else if (status != eslOK)   ESL_FAIL(status, errbuf, "Writing alignment file failed with error %d\n", status);
      }

      /* if nec, print #=GC RF annotation as a 1/0 mask (single line) to a file */
      if(esl_opt_GetString(go, "--omask") != NULL)
	{
	  if ((omaskfp = fopen(esl_opt_GetString(go, "--omask"), "w")) == NULL) 
	    ESL_FAIL(eslFAIL, errbuf, "Failed to open --omask output file %s\n", esl_opt_GetString(go, "--omask"));
	  if((status = output_rf_as_mask(omaskfp, errbuf, msa)) != eslOK) goto ERROR;
	}
      esl_msa_Destroy(msa);
      if(othermsa != NULL) esl_msa_Destroy(othermsa);
    }

  /* If an msa read failed, we drop out to here with an informative status code. 
   */
  if      (status == eslEFORMAT) 
    esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", 
	      afp->linenumber, afp->fname, afp->errbuf, afp->buf);	
  else if (status != eslEOF)
    esl_fatal("Alignment file read failed with error code %d\n", status);
  else if (nali   == 0)
    esl_fatal("No alignments found in file %s\n", alifile);

  if(esl_opt_GetString(go, "--omask") != NULL) fclose(omaskfp);
  if(esl_opt_GetString(go, "--kmask") != NULL) fclose(kmaskfp);
  
  /* Cleanup, normal return
   */
  if(otherafp != NULL) esl_msafile_Close(otherafp);
  if((esl_opt_GetString(go, "--morph") != NULL) && othermsa != NULL) esl_msa_Destroy(othermsa);

  if(esl_opt_GetString(go, "-o") != NULL) { fclose(ofp); }
  esl_msafile_Close(afp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  
  return 0;
  
 ERROR:
  if(afp != NULL) esl_msafile_Close(afp);
  if(go  != NULL) esl_getopts_Destroy(go);
  if(msa != NULL) esl_msa_Destroy(msa);

  esl_fatal(errbuf);
  return 1; /* never reached */
}


/* keep_or_remove_rf_gaps
 *                   
 * Given an MSA with #=GC RF markup, either remove or keep
 * all non-gap RF columns.
 */
static int
keep_or_remove_rf_gaps(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, int keep_flag, int remove_flag)
{
  int     status;
  int    *useme;
  int64_t apos;

  /* contract check */
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment.");
  if(keep_flag == TRUE && remove_flag == TRUE) ESL_XFAIL(eslEINVAL, errbuf, "in keep_or_remove_rf_gaps, keep_flag and remove_flag both TRUE.");
  if(keep_flag == FALSE && remove_flag == FALSE) ESL_XFAIL(eslEINVAL, errbuf, "in keep_or_remove_rf_gaps, keep_flag and remove_flag both FALSE.");

  ESL_ALLOC(useme, sizeof(int) * msa->alen);
  if(keep_flag)
  {
    for(apos = 0; apos < msa->alen; apos++)    
      useme[apos] = (esl_abc_CIsGap(msa->abc, msa->rf[apos]) ? FALSE : TRUE);
  }
  else if(remove_flag)
  {
    for(apos = 0; apos < msa->alen; apos++)    
      useme[apos] = (esl_abc_CIsGap(msa->abc, msa->rf[apos]) ? TRUE : FALSE);
  }
  else ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "In keep_or_remove_rf_gaps, but neither -r nor -k enabled.");
  if((status = esl_msa_ColumnSubset(msa, errbuf, useme)) != eslOK) return status;
  free(useme);
  return eslOK;

 ERROR:
  if(useme != NULL) free(useme);
  return eslEMEM;
}

/* keep_contiguous_column_block
 *                   
 * Keep only columns in range --start-all..--end-all, or --start-rf..--end-rf 
 */
static int
keep_contiguous_column_block(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa)
{
  int     status;
  int    *useme;
  int64_t apos;
  int     rf_mode;
  int     all_mode;
  int    *c2a_map = NULL;       /* msa map of consensus columns (non-gap RF residues) to alignment columns */
  int    clen;
  int    astart, aend;

  rf_mode  = ((!esl_opt_IsDefault(go, "--start-rf"))  && (!esl_opt_IsDefault(go, "--end-rf"))) ? TRUE : FALSE;
  all_mode = ((!esl_opt_IsDefault(go, "--start-all")) && (!esl_opt_IsDefault(go, "--end-all"))) ? TRUE : FALSE;
  if((!rf_mode) && (!all_mode)) ESL_XFAIL(eslEINVAL, errbuf, "Entered keep_contiguous_column_block, but neither (--start-rf & --end-rf) nor (--start-all & --end-all) combination invoked.");
  
  /* contract check */
  if(rf_mode && msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "--start-rf and --end-rf required #=GC RF markup in alignment, but none exists.");

  if(rf_mode) { 
    if((status = map_cpos_to_apos(msa, &c2a_map, &clen))   != eslOK) goto ERROR;
    if(esl_opt_GetInteger(go, "--start-rf") < 1)    ESL_XFAIL(eslEINVAL, errbuf, "<n> from --start-rf must be > 1.");
    if(esl_opt_GetInteger(go, "--end-rf")   > clen) ESL_XFAIL(eslEINVAL, errbuf, "<n> from --end-rf must be <= %d (which is the number of non-gap RF columns in the MSA).", clen);
    astart = c2a_map[esl_opt_GetInteger(go, "--start-rf")];
    aend   = c2a_map[esl_opt_GetInteger(go, "--end-rf")];
    if(astart > aend) ESL_XFAIL(eslEINVAL, errbuf, "<n> from --start-rf <n> must be lower than <n> from --end-rf.");
  }
  else { 
    if(esl_opt_GetInteger(go, "--start-all") < 1)         ESL_XFAIL(eslEINVAL, errbuf, "<n> from --start-all must be > 1.");
    if(esl_opt_GetInteger(go, "--end-all")   > msa->alen) ESL_XFAIL(eslEINVAL, errbuf, "<n> from --end-all must be <= %" PRId64 " (which is the number of columns in the MSA).", msa->alen);
    astart = esl_opt_GetInteger(go, "--start-all");
    aend   = esl_opt_GetInteger(go, "--end-all");
    if(astart > aend) ESL_XFAIL(eslEINVAL, errbuf, "<n> from --start-all <n> must be lower than <n> from --end-all.");
  }

  ESL_ALLOC(useme, sizeof(int) * msa->alen);
  esl_vec_ISet(useme, msa->alen, FALSE);
  for(apos = astart-1; apos < aend; apos++) useme[apos] = TRUE;
  if((status = esl_msa_ColumnSubset(msa, errbuf, useme)) != eslOK) return status;
  free(useme);
  if(c2a_map != NULL) free(c2a_map);
  return eslOK;

 ERROR:
  if(useme != NULL) free(useme);
  return eslEMEM;
}

/* write_rf_gapthresh
 *                   
 * Given an MSA write/rewrite RF based on fraction
 * of gaps in each column. If fraction < gapthresh RF is an 'x',
 * otherwise it's a '.' (gap).
 */
static int
write_rf_gapthresh(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa)
{
  int      status;
  int64_t  apos;
  int64_t  gaps;
  int      i;
  double   gapthresh;
  int      nrf = 0;

  if(msa->rf == NULL) { 
    ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen+1));
    for (apos = 1; apos <= msa->alen; apos++) msa->rf[(apos-1)] = '.';
  }

  gapthresh = esl_opt_GetReal(go, "--gapthresh");
  for (apos = 1; apos <= msa->alen; apos++)
    {
      for (gaps = 0, i = 0; i < msa->nseq; i++)
	if (esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) gaps++;
      if((double) gaps / (double) msa->nseq < gapthresh) { /* column passes gap threshold */
	nrf++;
	if(esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) msa->rf[(apos-1)] = 'x';
	/* else, leave it alone! */
      }
      else { /* column fails the gap threshold */
	msa->rf[(apos-1)] = '.';
      }
    }
  msa->rf[msa->alen] = '\0';

  if(esl_opt_GetBoolean(go, "--verbose")) printf("gapthresh %.3f %d of %d pass", gapthresh, nrf, (int) msa->alen);
  
  return eslOK;
 ERROR:
  return status;
}

/* write_rf_given_alen
 *                   
 * Given an MSA and a char string of 1s and 0s (a lanemask) of length
 * msa->alen, write/rewrite  RF positions as  'x' (non-gap) for 1, '.' (gap) for 0.
 * If RF already exists, do not modify non-gap RF columns if they are within mask ('1').
 */
static int
write_rf_given_alen(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, char *amask, int amask_len)
{
  int      status;
  int64_t  apos;

  /* contract check, rfgiven_mask must be exact length of msa */
  if(amask == NULL) ESL_FAIL(eslEINVAL, errbuf, "--mask-all mask is NULL in write_rf_given, this shouldn't happen.\n");
  if(amask_len != (int) strlen(amask)) { ESL_FAIL(eslEINVAL, errbuf, "write_rf_given_alen(), passed in mask len (%d) is not equal to actual mask length (%d)\n", amask_len, (int) strlen(amask)); }
  if(amask_len != msa->alen) 
    ESL_FAIL(eslEINVAL, errbuf, "--mask-all mask length: %d is not equal to the MSA length (%" PRId64 ")\n", 
	     amask_len, msa->alen); 
  if(msa->rf == NULL) { 
    ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen+1));
    for (apos = 1; apos <= msa->alen; apos++) msa->rf[(apos-1)] = '.';
  }

  for (apos = 1; apos <= msa->alen; apos++) {
    if     (amask[(apos-1)] == '0') msa->rf[(apos-1)] = '.';
    else if(amask[(apos-1)] == '1') { 
      if(esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) msa->rf[(apos-1)] = 'x'; /* else, leave it alone */
    }
    else    ESL_FAIL(eslEINVAL, errbuf, "--mask-all mask char number %" PRId64 " is not a 1 nor a 0, but a %c\n", apos, amask[(apos-1)]);
  }

  msa->rf[msa->alen] = '\0';
  return eslOK;
 ERROR:
  return status;
}

/* write_rf_given_rflen
 *
 * Given an MSA and a char string of 1s and 0s (a lanemask) that is
 * the same length as the non-gap RF annotation in msa, rewrite msa
 * RF based as 'x' (non-gap) for 1, '.' (gap) for 0. 1s indicate which
 * non-gap RF columns to keep as 'x', and 0s indicate which non-gap
 * RF columns to make gaps '.'.
 * If RF already exists, do not modify non-gap RF columns if they are 
 * within mask ('1').
 */
static int
write_rf_given_rflen(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, char *rfmask, int rfmask_len)
{
  int64_t  apos, cpos;

  /* contract check, mask must be exact length of msa */
  if(rfmask  == NULL) ESL_FAIL(eslEINVAL, errbuf, "--mask-rf mask is NULL in write_rf_given, this shouldn't happen.\n");
  if(msa->rf == NULL) ESL_FAIL(eslEINVAL, errbuf, "--mask-rf mask requires RF annotation in MSA (try -g)\n");
  if(rfmask_len != (int) strlen(rfmask)) { ESL_FAIL(eslEINVAL, errbuf, "write_rf_given_rflen(), passed in mask len (%d) is not equal to actual mask length (%d).\n", rfmask_len, (int) strlen(rfmask)); }

  cpos = 0;
  for (apos = 1; apos <= msa->alen; apos++) {
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) {
      cpos++;
      if     (rfmask[(cpos-1)] == '0') msa->rf[(apos-1)] = '.';
      else if(rfmask[(cpos-1)] == '1') { 
	if(esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) msa->rf[(apos-1)] = 'x'; /* else, leave it alone */
      }
    }
    else msa->rf[(apos-1)] = '.'; 
  }
  if(cpos != rfmask_len) { ESL_FAIL(eslEINVAL, errbuf, "write_rf_given_rflen(), RF non-gap length (consensus length) (%" PRId64 ") is not equal to mask length (%d)\n", cpos, rfmask_len); }

  msa->rf[msa->alen] = '\0';
  return eslOK;
}


/* write_rf_given_useme
 *
 * Given an MSA and a integer array <useme> of size msa->alen, set
 * msa->rf column [0..alen-1] i as 'x' if useme[i] == TRUE, and as
 * '.' if useme[i] == FALSE.
 */
static int
write_rf_given_useme(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, int *useme)
{
  int     status;
  int64_t apos;

  if(msa->rf == NULL) { 
    ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen+1));
    for (apos = 1; apos <= msa->alen; apos++) msa->rf[(apos-1)] = '.';
  }

  for (apos = 0; apos < msa->alen; apos++) { 
    if(useme[apos]) { 
      if(esl_abc_CIsGap(msa->abc, msa->rf[apos])) msa->rf[apos] = 'x'; /* else leave it alone */
    }
    else msa->rf[apos] = '.';
  }
  msa->rf[msa->alen] = '\0';
  return eslOK;
 ERROR:
  return status;
}

/* individualize_consensus
 *                   
 * Given an MSA with a consensus structure impose it to create
 * individual secondary structures. Simple rule, for consensus
 * bp i,j if seq positions i and j are both non-gaps seq i,j are 
 * paired, if >= 1 is a gap, they're not paired.
 */
static int
individualize_consensus(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa)
{
  int     status;
  int64_t apos;
  int  i;
  int *cct;		/* 0..alen-1 base pair partners array for consensus        */
  int *ct;		/* 0..alen-1 base pair partners array for current sequence */
  char *ss;             /* individual secondary structure we've built              */
  char *ss_cons_nopseudo; /* no-pseudoknot version of consensus structure */

  if(msa->ss_cons == NULL)                                ESL_FAIL(eslEINVAL, errbuf, "-i requires MSA to have consensus structure annotation.\n");
  if(! (msa->flags & eslMSA_DIGITAL))                     ESL_FAIL(eslEINVAL, errbuf, "individualize_consensus() MSA is not digitized.\n");
    
  ESL_ALLOC(cct, sizeof(int)  * (msa->alen+1));
  ESL_ALLOC(ct,  sizeof(int)  * (msa->alen+1));
  ESL_ALLOC(ss,  sizeof(char) * (msa->alen+1));
  ESL_ALLOC(ss_cons_nopseudo, sizeof(char) * (msa->alen+1));

  esl_wuss_nopseudo(msa->ss_cons, ss_cons_nopseudo);
  if (esl_wuss2ct(ss_cons_nopseudo, msa->alen, cct) != eslOK)     ESL_FAIL(status, errbuf, "Consensus structure string is inconsistent.");

  /* go through each position of each sequence, 
     if it's a gap and it is part of a base pair, remove that base pair */
  for (i = 0; i < msa->nseq; i++)
    {
      esl_vec_ICopy(cct, (msa->alen+1), ct);
      for (apos = 1; apos <= msa->alen; apos++)
	if (esl_abc_XIsGap(msa->abc, msa->ax[i][apos]))
	  { 
	    if (ct[apos] != 0)  ct[ct[apos]] = 0;
	    ct[apos] = 0;
	  }
      /* convert to WUSS SS string and append to MSA */
      if (esl_ct2wuss(ct, msa->alen, ss) != eslOK) ESL_FAIL(status, errbuf, "Consensus structure string had pseudoknots, we can't handle this yet.");
      esl_msa_AppendGR(msa, "SS", i, ss);
    }
  free(cct);
  free(ct);
  free(ss);
  free(ss_cons_nopseudo);
  return eslOK;
 ERROR:
  return status;
}

/* merge_msa
 *                   
 * Use the RF line as denoting consensus columns to merge
 * msa1 and msa2. msa1 is rewritten with merged msa. 
 * Important: msa1 will only contain sequence data from msa2.
 */
static int
merge_msa(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa1, ESL_MSA *msa2, ESL_MSA **ret_merged_msa)
{
  int status;
  int *agaps1 = NULL;
  int *agaps2 = NULL;
  int *c2a_map1 = NULL;       /* msa1 map of consensus columns (non-gap RF residues) to alignment columns */
  int *c2a_map2 = NULL;       /* msa2 map of consensus columns (non-gap RF residues) to alignment columns */


  int *new_c2a_map1 = NULL;   /* merged msa1 map of consensus columns (non-gap RF residues) to alignment columns */
  int *new_c2a_map2 = NULL;   /* merged msa2 map of consensus columns (non-gap RF residues) to alignment columns */
  int *aadd1 = NULL;          /* [1..apos..msa1->alen] number of columns to add after column apos to merge msa1 */
  int *aadd2 = NULL;          /* [1..apos..msa2->alen] number of columns to add after column apos to merge msa2 */
  int apos1, apos2;           /* counters over alignment positions of msa1, msa2 */
  int cpos = 0;               /* counter over consensus positions */
  int clen, clen2, new_clen1, new_clen2; /* consensus lengths */
  int tmp_ngaps;              /* temp var for number of gaps */
  int cur_apos1, nxt_apos1, cur_apos2, nxt_apos2; /* impt alignment positions */
  int astart2;                /* impt alignment positions */
  int ngaps1, ngaps2;         /* [1..apos..msa->alen] number of gaps in column apos of msa1, msa2 */
  int *msa2_cols_to_keep   = NULL; /* temp array for picking columns to keep in msa2 */
  int nadd1, nadd2;           /* temp vars, number of columns to keep, add */
  int radd = 0;               /* number of residues in msa2 columns corresponding to all 100% gap columns added to msa1 */
  int i, ip;                  /* sequence index counters */
  int x;                      /* general counter */
  int orig_msa1_nseq;         /* number of sequences in original msa1 */

  ESL_MSA *new_msa1;
  ESL_MSA *new_msa2;

  /* contract check */
  if(msa1->abc->type  != msa2->abc->type) ESL_XFAIL(eslEINVAL, errbuf, "With --merge both MSAs must have same alphabet.");
  if(msa1->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "With --merge both MSAs must have RF annotation.");
  if(msa2->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "With --merge both MSAs must have RF annotation.");
  
  /* Determine number of gaps in each column of msa1 and msa2 */
  if((status = get_gaps_per_column(msa1, &agaps1)) != eslOK) goto ERROR;
  if((status = get_gaps_per_column(msa2, &agaps2)) != eslOK) goto ERROR;

  /* Map consensus columns to alignment positions */
  if((status = map_cpos_to_apos(msa1, &c2a_map1, &clen))   != eslOK) goto ERROR;
  if((status = map_cpos_to_apos(msa2, &c2a_map2, &clen2))  != eslOK) goto ERROR;
  if(clen != clen2)
    ESL_XFAIL(eslEINVAL, errbuf, "With --merge both MSAs must have same consensus (non-gap RF) length.");
  
  /* Fill 'aadd1' and 'aadd2' arrays, these are the number columns of 100% gaps
   * we have to add after each 'msa1' and 'msa2' column respectively to make
   * the alignments the same size with identical non-gap RF lines.
   * Identical non-gap RF lines means for each aligned position i:
   *   isgap(msa1->rf[i]) == isgap(msa2->rf[i])
   */
  ESL_ALLOC(aadd1,  sizeof(int) * (msa1->alen+1));
  ESL_ALLOC(aadd2,  sizeof(int) * (msa2->alen+1));
  esl_vec_ISet(aadd1,  msa1->alen+1, 0);
  esl_vec_ISet(aadd2,  msa1->alen+1, 0);
  for(cpos = 0; cpos <= clen; cpos++)
    {
      if(cpos > 0) {
	cur_apos1 = c2a_map1[cpos];
	cur_apos2 = c2a_map2[cpos];
      }
      else cur_apos1 = cur_apos2 = 1;
      if(cpos < clen) {
	nxt_apos1 = c2a_map1[(cpos+1)];
	nxt_apos2 = c2a_map2[(cpos+1)];
      }
      else {
	nxt_apos1 = msa1->alen + 1;
	nxt_apos2 = msa2->alen + 1;
      }
      ngaps1 = nxt_apos1 - cur_apos1 - 1;
      ngaps2 = nxt_apos2 - cur_apos2 - 1;

      if(esl_opt_GetBoolean(go, "--verbose")) printf("%4d: ", cpos); 
      if(ngaps1 == ngaps2) /* we don't have to add any columns to either msa (okay if 0) */
	{
	  /* do nothing */
	  if(esl_opt_GetBoolean(go, "--verbose")) printf("\n");      
	}
      else if(ngaps1 <  ngaps2) /* we need to add some new 100% gap columns to msa1 */
	{ 
	  if(esl_opt_GetBoolean(go, "--verbose")) printf("\tmsa1 add     %4d all gap columns\n", (ngaps2-ngaps1)); 
	  nadd1 = ngaps2 - ngaps1;
	  /* determine where to put the gaps */
	  if(nxt_apos1 == (cur_apos1 + 1)) /* no choice, we have to put 100% gaps after cur_apos1 */
	    {
	      if(cpos == 0) aadd1[0] += nadd1;
	      else aadd1[c2a_map1[cpos]] += nadd1;
	    }
	  else
	    {
	      if(cpos == 0) 
		{ apos1 = astart2 = 0; }
	      else 
		{ 
		  apos1 = c2a_map1[cpos] + 1;
		  astart2 = cur_apos2+1; 
		}
	      tmp_ngaps = pick_gappiest_columns(agaps2, astart2, nxt_apos2-1, nadd1, &(msa2_cols_to_keep));
	      radd += (msa2->nseq * nadd1) - tmp_ngaps;
	      if(esl_opt_GetBoolean(go, "--verbose")) printf("\t\tresidues added: %d (%d)\n", ((msa2->nseq * nadd1) - tmp_ngaps), radd);
	      for(apos2 = astart2; apos2 < nxt_apos2; apos2++) 
		{
		  if(msa2_cols_to_keep[(apos2 - astart2)] == TRUE) aadd1[apos1]++;
		  else apos1++;
		}
	      if(apos1 != nxt_apos1) 
		esl_fatal("Coding error!");
	      free(msa2_cols_to_keep);
	    }
	}
      else if(ngaps1 >  ngaps2) /* we need to add some new 100% gap columns to msa 2 */
	{ 
	  if(esl_opt_GetBoolean(go, "--verbose")) printf("\tmsa2 add     %4d all gap columns\n", (ngaps1 - ngaps2));
	  nadd2 = ngaps1 - ngaps2;
	  /* determine where to put the gaps */
	  if(nxt_apos2 == (cur_apos2 + 1)) /* no choice, we have to put 100% gaps after cur_apos2 */
	    {
	      if(cpos == 0) aadd2[0] += nadd2;
	      else aadd2[c2a_map2[cpos]] += nadd2;
	    }
	  /*else
	    {
	      if(cpos == 0) 
		{ apos2 = astart1 = 0; }
	      else 
		{ 
		  apos2 = c2a_map2[cpos] + 1;
		  astart1 = cur_apos1+1; 
		}
	      tmp_ngaps = pick_gappiest_columns(agaps1, astart1, nxt_apos1-1, nadd, &(msa1_cols_to_keep));
	      radd += (msa2->nseq * nadd) - tmp_ngaps;
	      if(esl_opt_GetBoolean(go, "--verbose")) printf("\t\tresidues added: %d (%d)\n", ((msa2->nseq * nadd) - tmp_ngaps), radd);
	      for(apos2 = astart2; apos2 < nxt_apos2; apos2++) 
		{
		  if(msa2_cols_to_keep[(apos2 - astart2)] == TRUE) aadd1[apos1]++;
		  else apos1++;
		}
	      if(apos1 != nxt_apos1) 
		esl_fatal("Coding error!");
	      free(msa2_cols_to_keep);
	      }*/


	  /*if(cpos == 0) astart1 = 0;
	  else astart1 = cur_apos1+1;
	  if(ngaps2 == 0)
	    {
	      for(apos1 = astart1; apos1 < nxt_apos1; apos1++) akeep[apos1] = FALSE;
	      }*/
	  //else /* determine if it's likely flush left, or flush right first */
	  /*{
	      if(is_flush_left(agaps1, astart1, nxt_apos1-1))
		{
		  for(apos1 = astart1;         apos1 < (astart1 + nkeep); apos1++) 
		    akeep[apos1] = TRUE;
		  for(apos1 = (astart1 + nkeep); apos1 <  nxt_apos1;              apos1++)
		    akeep[apos1] = FALSE;
		}		  
	      else if(is_flush_right(agaps1, astart1, nxt_apos1-1))
		{
		  for(apos1 = astart1;             apos1 < (nxt_apos1 - nkeep); apos1++) 
		    akeep[apos1] = FALSE;
		  for(apos1 = (nxt_apos1 - nkeep); apos1 <  nxt_apos1;              apos1++)
		    akeep[apos1] = TRUE;
		    }*/
	  //else /* not flush left or flush right, pick least gappy columns to keep */
	  /*{
		  pick_gappiest_columns(agaps1, astart1, (nxt_apos1-1), (ngaps1 - nkeep), &(msa1_cols_to_remove));
		  for(apos1 = astart1; apos1 < nxt_apos1; apos1++) 
		    akeep[apos1] = (msa1_cols_to_remove[apos1 - astart1] == TRUE) ? FALSE : TRUE; 
		  free(msa1_cols_to_remove);
		}		
		}*/
	}
    }

  nadd1 = 0;
  if(esl_opt_GetBoolean(go, "--verbose")) { printf("Printing number of all gap columns to add after each msa1 alignment column:\n"); }
  for(apos1 = 1; apos1 <= msa1->alen; apos1++)
    {
      nadd1 += aadd1[apos1];
      if(esl_opt_GetBoolean(go, "--verbose")) { printf("%5d %5d\n", apos1, aadd1[apos1]); }
    }
  nadd1 += aadd1[0];
  if(esl_opt_GetBoolean(go, "--verbose")) printf("Adding  %d columns to msa 1\n", nadd1);

  nadd2 = 0;
  if(esl_opt_GetBoolean(go, "--verbose")) { printf("Printing number of all gap columns to add after each msa2 alignment column:\n"); }
  for(apos2 = 1; apos2 <= msa2->alen; apos2++)
    {
      nadd2 += aadd2[apos2];
      if(esl_opt_GetBoolean(go, "--verbose")) { printf("%5d %5d\n", apos2, aadd2[apos2]); }
    }
  nadd2 += aadd2[0];
  if(esl_opt_GetBoolean(go, "--verbose")) printf("Adding  %d columns to msa 2\n", nadd2);

  /* add the 100% gap columns to msa1 and msa2 */
  status = add_gap_columns_to_msa(errbuf, msa1, aadd1, &new_msa1, TRUE);
  status = add_gap_columns_to_msa(errbuf, msa2, aadd2, &new_msa2, TRUE);

  /* Make new_c2a_map1 and new_c2a_map2, they should be identical */
  if((status = map_cpos_to_apos(new_msa1, &new_c2a_map1, &new_clen1))  != eslOK) goto ERROR;
  if((status = map_cpos_to_apos(new_msa2, &new_c2a_map2, &new_clen2))  != eslOK) goto ERROR;
  if(new_clen1 != new_clen2) 
    ESL_XFAIL(eslEINVAL, errbuf, "Coding error, during alignment merge, after adding gaps, MSA lengths differ.");

  if(esl_opt_GetBoolean(go, "--verbose")) printf("printing final test\n\n");
  for(cpos = 1; cpos <= clen; cpos++) 
    {
      if(new_c2a_map1[cpos] != new_c2a_map2[cpos]) 
	esl_fatal("Coding error. Alignments to merge do not have same consensus position map\n");
      if(esl_opt_GetBoolean(go, "--verbose")) printf("%4d %4d %4d\n", cpos, new_c2a_map1[cpos], new_c2a_map2[cpos]);
    }

  /* merge msa2 into msa1 */


  /* first make sure all the info that should be the same is the same */
  if(new_msa1->alen  != new_msa2->alen)  esl_fatal("Coding error. Alignments to merge do not have same lengths.\n");
  if(new_msa1->flags != new_msa2->flags) esl_fatal("Alignments to merge do not have flags (this *could* be worked around, implement it if you want).\n");
  if(new_msa1->abc->type != new_msa2->abc->type) esl_fatal("Alignments to merge do not have same alphabet.\n");
  for(x = 0; x < eslMSA_NCUTS; x++) 
    {
      if     ( new_msa1->cutset[x] && !new_msa2->cutset[x]) esl_fatal("Alignments to merge do not have same cutoff info.\n");
      else if(!new_msa1->cutset[x] &&  new_msa2->cutset[x]) esl_fatal("Alignments to merge do not have same cutoff info.\n");
      else if( new_msa1->cutset[x] &&  new_msa2->cutset[x])
	if(fabs(new_msa1->cutoff[x] - new_msa2->cutoff[x]) > 0.0001)
	  esl_fatal("Alignments to merge do not have same cutoff info.\n");
    }

  /* now merge new_msa1 and new_msa2, by expanding new_msa1, and swapping ptrs to data in new_msa2, 
   * new_msa1 becomes merged alignment
   * new_msa2 becomes pathetic shell of an alignment
   *
   * to expand a MSA, the alen must be 0 (flag for esl_msa_Expand()) I 
   * reset it to -1 here and then back again. This may be ill advised 
   */
  new_msa1->alen = -1;
  while(new_msa1->sqalloc < (new_msa1->nseq + new_msa2->nseq)) 
    esl_msa_Expand(new_msa1);
  new_msa1->alen = new_msa2->alen;
  orig_msa1_nseq = new_msa1->nseq;
  
  if((new_msa1->ss_cons == NULL && new_msa2->ss_cons != NULL) ||
     (new_msa1->ss_cons != NULL && new_msa2->ss_cons == NULL) ||
     ((new_msa1->ss_cons != NULL && new_msa2->ss_cons != NULL) && 
      (strcmp(new_msa1->ss_cons, new_msa2->ss_cons)   != 0))) esl_fatal("Alignments to merge do not have same consensus structure.\n");
  if((new_msa1->sa_cons == NULL && new_msa2->sa_cons != NULL) ||
     (new_msa1->sa_cons != NULL && new_msa2->sa_cons == NULL) ||
     ((new_msa1->sa_cons != NULL && new_msa2->sa_cons != NULL) && 
      (strcmp(new_msa1->sa_cons, new_msa2->sa_cons)   != 0))) esl_fatal("Alignments to merge do not have same consensus structure.\n");
  if((new_msa1->aseq == NULL && new_msa2->aseq != NULL) ||
     (new_msa1->aseq != NULL && new_msa2->aseq == NULL))  esl_fatal("Alignments to merge aseqs null/non-null mismatch.\n");
#ifdef eslAUGMENT_ALPHABET
  if((new_msa1->ax == NULL && new_msa2->ax != NULL) ||
     (new_msa1->ax != NULL && new_msa2->ax == NULL))  esl_fatal("Alignments to merge ax null/non-null mismatch.\n");
#endif /*eslAUGMENT_ALPHABET*/
  if((new_msa1->sqacc == NULL && new_msa2->sqacc != NULL) ||
     (new_msa1->sqacc != NULL && new_msa2->sqacc == NULL))  esl_fatal("Alignments to merge sqacc null/non-null mismatch.\n");
  if((new_msa1->sqdesc == NULL && new_msa2->sqdesc != NULL) ||
     (new_msa1->sqdesc != NULL && new_msa2->sqdesc == NULL))  esl_fatal("Alignments to merge sqdesc null/non-null mismatch.\n");
  if((new_msa1->ss == NULL && new_msa2->ss != NULL) ||
     (new_msa1->ss != NULL && new_msa2->ss == NULL))  esl_fatal("Alignments to merge ss null/non-null mismatch.\n");
  if((new_msa1->sa == NULL && new_msa2->sa != NULL) ||
     (new_msa1->sa != NULL && new_msa2->sa == NULL))  esl_fatal("Alignments to merge sa null/non-null mismatch.\n");

     /* rf lines were already indirectly checked, they must be equal */

  for(i = orig_msa1_nseq; i < (orig_msa1_nseq + new_msa2->nseq); i++)
    {
      ip = i - orig_msa1_nseq;
      if(new_msa1->aseq != NULL) new_msa1->aseq[i]   = new_msa2->aseq[ip];
#ifdef eslAUGMENT_ALPHABET
      if(new_msa1->ax   != NULL) new_msa1->ax[i]     = new_msa2->ax[ip];
#endif /*eslAUGMENT_ALPHABET*/
      new_msa1->sqname[i] = new_msa2->sqname[ip];
      new_msa1->wgt[i]    = new_msa2->wgt[ip];
      new_msa1->nseq++;

      if(new_msa1->sqacc  != NULL) new_msa1->sqacc[i]  = new_msa2->sqacc[ip];
      if(new_msa1->sqdesc != NULL) new_msa1->sqdesc[i] = new_msa2->sqdesc[ip];
      if(new_msa1->ss     != NULL) new_msa1->ss[i]     = new_msa2->ss[ip];
      if(new_msa1->sa     != NULL) new_msa1->sa[i]     = new_msa2->sa[ip];

      /* new_msa1->name,desc,acc,au untouched (is this unwise?) */

      if(new_msa1->sqlen  != NULL) new_msa1->sqlen[i]  = new_msa2->sqlen[ip];

      if(new_msa1->sslen != NULL) new_msa1->sslen[i]  = new_msa2->sslen[ip];
      if(new_msa1->salen != NULL) new_msa1->salen[i]  = new_msa2->salen[ip];
      /* lastidx not touched, should be unimportant */
    }
  /* copy and free comments (no need to swap pointers thanks to convenient esl_msa_AddComment() function */
  for(x = 0; x < new_msa2->ncomment; x++) {
    esl_msa_AddComment(new_msa1, new_msa2->comment[x]);
    free(new_msa2->comment[x]);
  }
  /* copy and free GF markup */
  for(x = 0; x < new_msa2->ngf; x++) {
    esl_msa_AddGF(new_msa1, new_msa2->gf_tag[x], new_msa2->gf[x]);
    free(new_msa2->gf_tag[x]);
    free(new_msa2->gf[x]);
  }
  /* copy and free GS markup */
  for(x = 0; x < new_msa2->ngs; x++) {
    for(i = orig_msa1_nseq; i < (orig_msa1_nseq + new_msa2->nseq); i++)
      {
	ip = i - orig_msa1_nseq;
	esl_msa_AddGS(new_msa1, new_msa2->gs_tag[x], i, new_msa2->gs[x][ip]);
	free(new_msa2->gs[x][ip]);
      }
    free(new_msa2->gs_tag[x]);
  }
  /* don't touch GC (per column) annotation, is this unwise? */

  /* copy and free GR markup */
  for(x = 0; x < new_msa2->ngr; x++) {
    for(i = orig_msa1_nseq; i < (orig_msa1_nseq + new_msa2->nseq); i++)
      {
	ip = i - orig_msa1_nseq;
	esl_msa_AppendGR(new_msa1, new_msa2->gr_tag[x], i, new_msa2->gr[x][ip]);
	free(new_msa2->gr[x][ip]);
      }
    free(new_msa2->gr_tag[x]);
  }
  /* don't touch keyhashes or SSI offset, shouldn't be a problem since we're just
   * printing the alignment. 
   */

  *ret_merged_msa = new_msa1;
  free(new_msa2); /* the guts are still valid, being pointed to by new_msa1 */

  free(agaps1);
  free(agaps2);
  free(c2a_map1);
  free(c2a_map2);
  free(new_c2a_map1);
  free(aadd1);
  free(aadd2);
  return eslOK;

 ERROR:
  if(agaps1       != NULL) free(agaps1);
  if(agaps2       != NULL) free(agaps2);
  if(c2a_map1     != NULL) free(c2a_map1);
  if(c2a_map2     != NULL) free(c2a_map2);
  if(new_c2a_map1 != NULL) free(new_c2a_map1);
  if(aadd1        != NULL) free(aadd1);
  if(aadd2        != NULL) free(aadd2);
  return status;
}

/* morph_msa
 *                   
 * Use the RF line as denoting consensus columns to morph
 * msa1 into msa2's gap structure. This may require removing
 * some columns from msa1, and adding some 100% gap columns
 * to msa1.
 */
static int
morph_msa(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa1, ESL_MSA *msa2, ESL_MSA **new_msa1)
{
  int status;
  int *agaps1 = NULL;
  int *agaps2 = NULL;
  int *c2a_map1 = NULL;       /* msa1 map of consensus columns (non-gap RF residues) to alignment columns */
  int *c2a_map2 = NULL;       /* msa2 map of consensus columns (non-gap RF residues) to alignment columns */
  int *new_c2a_map1 = NULL;   /* morphed msa1 map of consensus columns (non-gap RF residues) to alignment columns */
  int *akeep = NULL;          /* [1..apos..msa1->alen] TRUE to keep column apos in morphed msa1, FALSE not to */
  int *aadd = NULL;           /* [1..apos..msa1->alen] number of columns to add after column apos to morphed msa1 */
  int apos1, apos2;           /* counters over alignment positions of msa1, msa2 */
  int cpos = 0;               /* counter over consensus positions */
  int clen, clen2, new_clen1; /* consensus lengths */
  int tmp_ngaps;              /* temp var for number of gaps */
  int cur_apos1, nxt_apos1, cur_apos2, nxt_apos2; /* impt alignment positions */
  int astart1, astart2;       /* impt alignment positions */
  int ngaps1, ngaps2;         /* [1..apos..msa->alen] number of gaps in column apos of msa1, msa2 */
  int *msa1_cols_to_remove = NULL; /* temp array for picking columns to remove from msa1 */
  int *msa2_cols_to_keep   = NULL; /* temp array for picking columns to keep in msa2 */
  int nkeep, nadd;            /* temp vars, number of columns to keep, add */
  int radd = 0;               /* number of residues in msa2 columns corresponding to all 100% gap columns added to msa1 */
  int delete = 0;             /* number of residues in msa1 we have to delete during morph */

  /* contract check */
  if(msa1->abc->type  != msa2->abc->type) ESL_XFAIL(eslEINVAL, errbuf, "With --morph both MSAs must have same alphabet.");
  if(msa1->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "With --morph both MSAs must have RF annotation.");
  if(msa2->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "With --morph both MSAs must have RF annotation.");
  
  /* Determine number of gaps in each column of msa1 and msa2 */
  if((status = get_gaps_per_column(msa1, &agaps1)) != eslOK) goto ERROR;
  if((status = get_gaps_per_column(msa2, &agaps2)) != eslOK) goto ERROR;

  /* Map consensus columns to alignment positions */
  if((status = map_cpos_to_apos(msa1, &c2a_map1, &clen))   != eslOK) goto ERROR;
  if((status = map_cpos_to_apos(msa2, &c2a_map2, &clen2))  != eslOK) goto ERROR;
  if(clen != clen2)
    ESL_XFAIL(eslEINVAL, errbuf, "With --morph both MSAs must have same consensus (non-gap RF) length.");
  
  /* Fill the 'akeep' array [1..msa1->alen] and 'aadd' array [0..msa1->alen] which 
   * tells us which columns in msa1 we'll keep, and how many columns of 100% gaps
   * to add after each msa1 column. 
   */
  ESL_ALLOC(akeep, sizeof(int) * (msa1->alen+1));
  ESL_ALLOC(aadd,  sizeof(int) * (msa1->alen+1));
  esl_vec_ISet(akeep, msa1->alen+1, FALSE);
  esl_vec_ISet(aadd,  msa1->alen+1, 0);
  for(cpos = 0; cpos <= clen; cpos++)
    {
      if(cpos > 0) {
	cur_apos1 = c2a_map1[cpos];
	cur_apos2 = c2a_map2[cpos];
      }
      else cur_apos1 = cur_apos2 = 1;
      if(cpos < clen) {
	nxt_apos1 = c2a_map1[(cpos+1)];
	nxt_apos2 = c2a_map2[(cpos+1)];
      }
      else {
	nxt_apos1 = msa1->alen + 1;
	nxt_apos2 = msa2->alen + 1;
      }
      akeep[cur_apos1] = TRUE; /* keep the consensus column */
      ngaps1 = nxt_apos1 - cur_apos1 - 1;
      ngaps2 = nxt_apos2 - cur_apos2 - 1;

      if(esl_opt_GetBoolean(go, "--verbose")) printf("%4d: ", cpos); 
      if(ngaps1 == ngaps2) /* keep all columns in between (okay if 0) */
	{
	  for(apos1 = cur_apos1+1; apos1 < nxt_apos1; apos1++) akeep[apos1] = TRUE; 
	  if(esl_opt_GetBoolean(go, "--verbose")) printf("\n");      
	}
      else if(ngaps1 <  ngaps2) /* we need to add some new 100% gap columns */
	{ 
	  if(esl_opt_GetBoolean(go, "--verbose")) printf("\tadd     %4d all gap columns\n", (ngaps2-ngaps1)); 
	  nadd = ngaps2 - ngaps1;
	  /* keep all the inserts we have in msa1 */
	  for(apos1 = cur_apos1+1; apos1 < nxt_apos1; apos1++) akeep[apos1] = TRUE;
	  if(nxt_apos1 == (cur_apos1 + 1)) /* no choice, we have to put 100% gaps after cur_apos1 */
	    {
	      if(cpos == 0) aadd[0] += nadd;
	      else aadd[c2a_map1[cpos]] += nadd;
	    }
	  else
	    {
	      if(cpos == 0) 
		{ apos1 = astart2 = 0; }
	      else 
		{ 
		  apos1 = c2a_map1[cpos] + 1;
		  astart2 = cur_apos2+1; 
		}
	      tmp_ngaps = pick_gappiest_columns(agaps2, astart2, nxt_apos2-1, nadd, &(msa2_cols_to_keep));
	      radd += (msa2->nseq * nadd) - tmp_ngaps;
	      if(esl_opt_GetBoolean(go, "--verbose")) printf("\t\tresidues added: %d (%d)\n", ((msa2->nseq * nadd) - tmp_ngaps), radd);
	      for(apos2 = astart2; apos2 < nxt_apos2; apos2++) 
		{
		  if(msa2_cols_to_keep[(apos2 - astart2)] == TRUE) aadd[apos1]++;
		  else apos1++;
		}
	      if(apos1 != nxt_apos1) 
		esl_fatal("Coding error 10.");
	      free(msa2_cols_to_keep);
	    }
	}
      else if(ngaps1 >  ngaps2) /* we need to delete some of our msa1 columns */
	{ 
	  if(esl_opt_GetBoolean(go, "--verbose")) printf("\tdelete  %4d/%4d    columns\n", (ngaps1 - ngaps2), (ngaps1));  
	  nkeep = ngaps2;
	  if(cpos == 0) astart1 = 0;
	  else astart1 = cur_apos1+1;
	  if(ngaps2 == 0)
	    {
	      for(apos1 = astart1; apos1 < nxt_apos1; apos1++) akeep[apos1] = FALSE;
	    }
	  else /* determine if it's likely flush left, or flush right first */
	    {
	      if(is_flush_left(agaps1, astart1, nxt_apos1-1))
		{
		  for(apos1 = astart1;         apos1 < (astart1 + nkeep); apos1++) 
		    akeep[apos1] = TRUE;
		  for(apos1 = (astart1 + nkeep); apos1 <  nxt_apos1;              apos1++)
		    akeep[apos1] = FALSE;
		}		  
	      else if(is_flush_right(agaps1, astart1, nxt_apos1-1))
		{
		  for(apos1 = astart1;             apos1 < (nxt_apos1 - nkeep); apos1++) 
		    akeep[apos1] = FALSE;
		  for(apos1 = (nxt_apos1 - nkeep); apos1 <  nxt_apos1;              apos1++)
		    akeep[apos1] = TRUE;
		}
	      else /* not flush left or flush right, pick least gappy columns to keep */
		{
		  pick_gappiest_columns(agaps1, astart1, (nxt_apos1-1), (ngaps1 - nkeep), &(msa1_cols_to_remove));
		  for(apos1 = astart1; apos1 < nxt_apos1; apos1++) 
		    akeep[apos1] = (msa1_cols_to_remove[apos1 - astart1] == TRUE) ? FALSE : TRUE; 
		  free(msa1_cols_to_remove);
		}		
	    }
	}
    }

  nadd = 0;
  nkeep = 0;
  if(esl_opt_GetBoolean(go, "--verbose")) { printf("Printing number of all gap columns to add after each msa1 alignment column:\n"); }
  for(apos1 = 1; apos1 <= msa1->alen; apos1++)
    {
      if(akeep[apos1]) nkeep++;
      else delete += (msa1->nseq - agaps1[apos1]);
      nadd += aadd[apos1];
      if(esl_opt_GetBoolean(go, "--verbose")) { printf("%5d %5d\n", apos1, aadd[apos1]); }
    }
  nadd += aadd[0];
  printf("\n\nKeeping %d columns, deleting %d residues.\n", nkeep, delete);
  printf("Adding  %d columns, which have %d total non-gaps in MSA2.\n", nadd, radd);

  /* Rewrite the msa->rf line so that we can call keep_or_remove_rf_gaps to remove
   * the columns we don't want. Then restore the rf line. Do this by making a new
   * #=GC ORIGRF markup line. This is a serious violation of easel conventions.
   */
  char *origrf;
  esl_strdup(msa1->rf, msa1->alen, &origrf);
  esl_msa_AppendGC(msa1, "ORIGRF", origrf);
  /* overwrite RF temporarily with 'x's for any column we're keeping, '.' for any we're losing */
  for(apos1 = 1; apos1 <= msa1->alen; apos1++)
    msa1->rf[(apos1-1)] = akeep[apos1] == FALSE ? '.' : 'x';
  
  /* add the 100% gap columns */
  status = add_gap_columns_to_msa(errbuf, msa1, aadd, new_msa1, FALSE);

  /* remove unwanted columns */
  keep_or_remove_rf_gaps(go, errbuf, *new_msa1, TRUE, FALSE); 

  /* restore RF line */
  free((*new_msa1)->rf);
  esl_strdup((*new_msa1)->gc[(*new_msa1)->ngc-1], (*new_msa1)->alen, &((*new_msa1)->rf));
  free(origrf);
  free((*new_msa1)->gc_tag[((*new_msa1)->ngc-1)]);
  free((*new_msa1)->gc[((*new_msa1)->ngc-1)]);
  (*new_msa1)->ngc--;

  /* Make new new_c2a_map1, it should be identical to c2a_map2. */
  if((status = map_cpos_to_apos((*new_msa1), &new_c2a_map1, &new_clen1))  != eslOK) goto ERROR;
  if(new_clen1 != clen) 
    ESL_XFAIL(eslEINVAL, errbuf, "With --morph both MSAs must have same consensus (non-gap RF) length.");

  if(esl_opt_GetBoolean(go, "--verbose")) printf("printing final test\n\n");
  for(cpos = 1; cpos <= clen; cpos++) 
    {
      if(c2a_map2[cpos] != new_c2a_map1[cpos]) 
	esl_fatal("Coding error. Morphed alignment does not have same consensus position map as %s\n", esl_opt_GetString(go, "--morph"));
      if(esl_opt_GetBoolean(go, "--verbose")) printf("%4d %4d %4d %4d\n", cpos, c2a_map2[cpos], new_c2a_map1[cpos], (c2a_map2[cpos] - new_c2a_map1[cpos]));
    }
  
  free(agaps1);
  free(agaps2);
  free(c2a_map1);
  free(c2a_map2);
  free(new_c2a_map1);
  free(akeep);
  free(aadd);
  return eslOK;

 ERROR:
  if(agaps1       != NULL) free(agaps1);
  if(agaps2       != NULL) free(agaps2);
  if(c2a_map1     != NULL) free(c2a_map1);
  if(c2a_map2     != NULL) free(c2a_map2);
  if(new_c2a_map1 != NULL) free(new_c2a_map1);
  if(akeep        != NULL) free(akeep);
  if(aadd         != NULL) free(aadd);
  return status;
}

/* add_gap_columns_to_msa
 *                   
 * Given an MSA and an array specifying a number
 * of all gap columns to add after each column,
 * add them. Reallocate all arrays as necessary.
 * if(do_treat_as_rf_gap) make new column a gap
 * in the RF line, else make it an 'x'.
 *
 * toadd is numbered 1..alen.
 */
static int
add_gap_columns_to_msa(char *errbuf, ESL_MSA *msa, int *toadd, ESL_MSA **ret_msa, int do_treat_as_rf_gap)
{
  int status;
  int i,j;
  int apos;
  int nnew = 0;
  ESL_ALPHABET *abc;
  char *newstr;
  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL))
    ESL_XFAIL(eslEINVAL, errbuf, "in add_gap_columns_to_msa(), msa must be digitized.");
  for(apos = 0; apos <= msa->alen; apos++)
    nnew += toadd[apos];

  /* Textize the alignment */
  abc = msa->abc;
  esl_msa_Textize(msa);

  ESL_MSA *newmsa;
  /*printf("msa->nseq: %d\n", msa->nseq);
    printf("msa->alen: %d\n", msa->alen);*/
  newmsa = esl_msa_Create(msa->nseq, (msa->alen+nnew));

  /* Copy and add gaps to all valid data that is [0..(alen-1)] or [1..alen] */ 
  if(msa->ss_cons != NULL) 
    {
      ESL_ALLOC(newmsa->ss_cons, sizeof(char) * (msa->alen+nnew+1));
      if((status = cp_and_add_gaps_to_aseq(newmsa->ss_cons, msa->ss_cons, msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
    }
  if(msa->sa_cons != NULL) 
    {
      ESL_ALLOC(newmsa->sa_cons, sizeof(char) * (msa->alen+nnew+1));
      if((status = cp_and_add_gaps_to_aseq(newmsa->sa_cons, msa->sa_cons, msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
    }
  if(msa->rf != NULL)
    {
      ESL_ALLOC(newmsa->rf, sizeof(char) * (msa->alen+nnew+1));
      if(do_treat_as_rf_gap)
	{
	  if((status = cp_and_add_gaps_to_aseq(newmsa->rf,      msa->rf,      msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
	}
      else if((status = cp_and_add_gaps_to_aseq(newmsa->rf,      msa->rf,      msa->alen, toadd, nnew, 'x') != eslOK)) goto ERROR;
    }

  if(msa->ss != NULL)
    {
      ESL_ALLOC(newmsa->ss, sizeof(char *) * msa->nseq);
      for(i = 0; i < msa->nseq; i++)
      {
	if(msa->ss[i] != NULL)
	  {
	    ESL_ALLOC(newmsa->ss[i], sizeof(char) * (msa->alen+nnew+1));
	    if((status = cp_and_add_gaps_to_aseq(newmsa->ss[i], msa->ss[i], msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
	  }
      }
    }

  if(msa->sa != NULL)
    {
      for(i = 0; i < msa->nseq; i++)
      {
	if(msa->sa[i] != NULL)
	  {
	    ESL_ALLOC(newmsa->sa[i], sizeof(char) * (msa->alen+nnew+1));
	    if((status = cp_and_add_gaps_to_aseq(newmsa->sa[i], msa->sa[i], msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
	  }
      }
    }  

  if(msa->ncomment > 0)
    {
      for(j = 0; j < msa->ncomment; j++)
	{
	  if(msa->comment[j] != NULL) 
	    esl_msa_AddComment(newmsa, msa->comment[j]);
	}
    }

  if(msa->ngf > 0)
    {
      for(i = 0; i < msa->ngf; i++)
	if(msa->gf[i] != NULL) 
	    esl_msa_AddGF(newmsa, msa->gf_tag[i], msa->gf[i]);
    }

  if(msa->ngs > 0)
    {
      for(j = 0; j < msa->ngs; j++)
	{
	  for(i = 0; i < msa->nseq; i++)
	    if(msa->gs[j][i] != NULL) 
	      esl_msa_AddGS(newmsa, msa->gs_tag[j], i, msa->gs[j][i]);
	}
    }

  if(msa->ngc > 0)
    {
      for(i = 0; i < msa->ngc; i++)
	{
	  if(msa->gc[i] != NULL) 
	    {
	      ESL_ALLOC(newstr, sizeof(char) * (msa->alen+nnew+1));
	      if((status = cp_and_add_gaps_to_aseq(newstr, msa->gc[i], msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
	      esl_msa_AppendGC(newmsa, msa->gc_tag[i], newstr);
	      free(newstr);
	    }
	}
    }

  if(msa->gr != NULL)
    {  
      for(j = 0; j < msa->ngr; j++)
	{
	  for(i = 0; i < msa->nseq; i++)
	    {
	      if(msa->gr[j][i] != NULL) 
		{
		  ESL_ALLOC(newstr, sizeof(char) * (msa->alen+nnew+1));
		  if((status = cp_and_add_gaps_to_aseq(newstr, msa->gr[j][i], msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
		  esl_msa_AppendGR(newmsa, msa->gr_tag[j], i, newstr);
		  free(newstr);
		}
	    }
	}
    }
    
  /* copy the aseqs, free as we go to save memory */
  for(i = 0; i < msa->nseq; i++)
    {
      esl_strdup(msa->sqname[i], -1, &(newmsa->sqname[i]));
      if((status = cp_and_add_gaps_to_aseq(newmsa->aseq[i], msa->aseq[i], msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
      free(msa->aseq[i]);
      msa->aseq[i] = NULL;
    }    
  newmsa->abc = abc;
  esl_msa_Digitize(newmsa->abc, newmsa);
  esl_msa_Destroy(msa);
  *ret_msa = newmsa;
      
  return eslOK;
  
 ERROR:
  return status;
}

/*cp_and_add_gaps_to_aseq
 *                   
 * Given an aligned [0..alen-1] original text string,
 * add toadd[apos-1] gaps after each residue. 
 * new_aseq must be already allocated. 
 *
 * toadd is numbered 1..alen.
 */
static int cp_and_add_gaps_to_aseq(char *new_aseq, char *orig_aseq, int alen, int *toadd, int nnew, char gapchar)
{
  int orig_apos = 0;
  int new_apos  = 0;
  int i;

  for(i = 0; i < toadd[0]; i++)
    new_aseq[new_apos++] = gapchar;
  for(orig_apos = 0; orig_apos < alen; orig_apos++)
    {
      new_aseq[new_apos++] = orig_aseq[orig_apos];
      for(i = 0; i < toadd[(orig_apos+1)]; i++)
	new_aseq[new_apos++] = gapchar;
    }
  new_aseq[new_apos] = '\0';
  return eslOK;
}

/* is_flush_left
 *                   
 * Given an array with number of gaps in each column
 * of an alignment, and an interval of columns astart..aend,
 * return TRUE if the residues in this interval appear to 
 * leftflushed inserts, FALSE otherwise
 */
static int is_flush_left(int *ngaps, int astart, int aend)
{
  if(astart == -1 || aend == -1) esl_fatal("is_flush_left invalid column positions.");
  
  int i;
  int gaps = ngaps[astart];
  for(i = astart+1; i <= aend; i++)
    {
      if(ngaps[i] < gaps) return FALSE;
      gaps = ngaps[i];
    }
  return TRUE;
}

/* is_flush_right
 *                   
 * Given an array with number of gaps in each column
 * of an alignment, and an interval of columns astart..aend,
 * return TRUE if the residues in this interval appear to 
 * rightflushed inserts, FALSE otherwise
 */
static int is_flush_right(int *ngaps, int astart, int aend)
{
  if(astart == -1 || aend == -1) esl_fatal("is_flush_right invalid column positions.");
  
  int i;
  int gaps = ngaps[astart];
  for(i = astart+1; i <= aend; i++)
    {
      if(ngaps[i] > gaps) return FALSE;
      gaps = ngaps[i];
    }
  return TRUE;
}

/* pick_gappiest_columns
 *                   
 * Given an array with number of gaps in each column
 * of an alignment, and an interval of columns astart..aend.
 * Pick the npick gappiest columns and store that info
 * in ret_cols_to_pick.
 *
 * Returns total number of gaps the npick columns picked.
 */
static int pick_gappiest_columns(int *ngaps, int astart, int aend, int npick, int **ret_cols_to_pick)
{
  if(astart == -1 || aend == -1) esl_fatal("pick_gappiest_columns invalid column positions.");
  if((aend-astart+1) < npick)    esl_fatal("pick_gappiest_columns number to pick (%d) exceeds number of possibilities (%d).", npick, (aend-astart+1));

  int status;
  int i,c;
  int *tmp_ngaps;
  int *cols_to_pick;
  int topick;
  int total_gaps = 0;

  ESL_ALLOC(tmp_ngaps,    sizeof(int) * (aend-astart+1));
  ESL_ALLOC(cols_to_pick, sizeof(int) * (aend-astart+1));

  esl_vec_ISet(cols_to_pick, (aend-astart+1), FALSE);
  for(i = astart; i <= aend; i++)
    tmp_ngaps[(i-astart)] = ngaps[astart];
  for(c = 0; c < npick; c++)
    {
      topick               = esl_vec_IArgMax(tmp_ngaps, (aend-astart+1));
      cols_to_pick[topick] = TRUE;
      total_gaps          += tmp_ngaps[topick];
      tmp_ngaps[topick]    = -1;
    }
  free(tmp_ngaps);
  *ret_cols_to_pick = cols_to_pick;
  return total_gaps;

 ERROR:
  esl_fatal("Memory allocation error.");
  return -1;
}

/* get_gaps_per_column 
 *                   
 * Given an MSA, determine the number of gaps per
 * column, and return a newly allocated array with this
 * into in *ret_ngaps. 
 */
static int get_gaps_per_column(ESL_MSA *msa, int **ret_ngaps)
{
  int status;
  int i, apos;
  int *ngaps = NULL;
  /* contract check */
  if(! msa->flags & eslMSA_DIGITAL) { status = eslEINVAL; goto ERROR; }

  ESL_ALLOC(ngaps, sizeof(int) * (msa->alen+1));
  esl_vec_ISet(ngaps, msa->alen+1, 0);
  for(i = 0; i < msa->nseq; i++) {
    for(apos = 1; apos <= msa->alen; apos++)
      ngaps[apos] += esl_abc_XIsGap(msa->abc, msa->ax[i][apos]);
  }
  *ret_ngaps = ngaps;
  return eslOK;

 ERROR:
  if(ngaps != NULL) free(ngaps);
  return status;
}

/* map_cpos_to_apos
 *                   
 * Given an MSA, determine the alignment position each
 * consensus position refers to. 
 */
static int map_cpos_to_apos(ESL_MSA *msa, int **ret_c2a_map, int *ret_clen)
{
  int status;
  int clen = 0;
  int *c2a_map = NULL;
  int cpos = 0;
  int apos = 0;
  /* contract check */
  if(msa->rf == NULL) { status = eslEINVAL; goto ERROR; }

  /* count consensus columns */
  for(apos = 1; apos <= msa->alen; apos++)
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) clen++;

  /* build map */
  ESL_ALLOC(c2a_map, sizeof(int) * (clen+1));
  c2a_map[0] = -1;
  for(apos = 1; apos <= msa->alen; apos++) 
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) c2a_map[++cpos] = apos;

  *ret_c2a_map = c2a_map;
  *ret_clen    = clen;
  return eslOK;

 ERROR:
  if(c2a_map != NULL) free(c2a_map);
  return status;
}


/* read_sqfile
 *                   
 * Read all seqs in a sequence file and return them. Originally
 * written for --trim option.
 */
static int read_sqfile(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc, int nseq, ESL_SQ ***ret_sq)
{
  int status;
  ESL_SQ **sq; 
  int i;
  
  /* get seqs from sqfile */
  ESL_ALLOC(sq, sizeof(ESL_SQ *) * (nseq + 1)); /* +1 for the last guy we allocate but don't use */
  i = 0;
  sq[i] = esl_sq_CreateDigital(abc);
  while ((status = esl_sqio_Read(sqfp, sq[i])) == eslOK) { 
    i++;
    if(i > nseq) esl_fatal("With --trim, sequence file must have same number seqs as in <msafile>\n"); 
    sq[i] = esl_sq_CreateDigital(abc);
  }
  if (i != nseq) esl_fatal("With --trim, sequence file must have same number seqs as in <msafile>\n"); 
  /* status should be eslEOF on normal end; if it isn't, deal w/ error */
  esl_sq_Destroy(sq[i]); /* destroy final allocated but unused seq */
    if (status == eslEFORMAT)
      esl_fatal("\
Sequence file parse error, line %d of file %s:\n\
%s\n", sqfp->linenumber, sqfp->filename, sqfp->errbuf);
    else if (status != eslEOF)
      esl_fatal("Sequence file %s read failed with error code %d\n",
		sqfp->filename, status);
  esl_sqfile_Close(sqfp);
  *ret_sq = sq;

  return eslOK;

 ERROR:
  esl_fatal("Memory allocation error.");
  return status; /* NEVERREACHED */
}


/* trim_msa
 *                   
 * Given an MSA and unaligned 'trimmed' versions (subsequences) of all seqs in that MSA, 
 * replace all chars that have been trimmed away (not in subsequences) with gaps in the MSA.
 */
static int trim_msa(ESL_MSA *msa, ESL_SQ **sq, char *errbuf)
{
  int status;
  int i;
  int apos, uapos;
  int astart,  aend;
  int uastart, uaend;
  char *offset;
  char *aseq;
  char *uaseq;
  char *uasubseq;
  int *a2ua_map;
  int *ua2a_map;
  int ualen;

  if(! (msa->flags & eslMSA_DIGITAL))
    ESL_XFAIL(eslEINVAL, errbuf, "in trim_msa(), msa must be digitized.");

  ESL_ALLOC(aseq,  sizeof(char) * (msa->alen+1));

  for(i = 0; i < msa->nseq; i++)
    {
      if (sq[i]->dsq == NULL) ESL_XFAIL(eslEINVAL, errbuf, "in trim_msa(), sq's must be digitized.");
      if (sq[i]->n   == 0)    ESL_XFAIL(eslEINVAL, errbuf, "in trim_msa(), sq[%d] is zero-length\n", i);

      ESL_ALLOC(a2ua_map, sizeof(int) * (msa->alen+1));
      esl_vec_ISet(a2ua_map, (msa->alen+1), -1);
      uapos = apos = 1;
      while(apos <= msa->alen)
	{
	  while(apos <= msa->alen && esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) apos++;
	  if(apos <= msa->alen) a2ua_map[apos] = uapos++;
	  apos++;
	}
      ualen = uapos;
      ESL_ALLOC(ua2a_map, sizeof(int) * (ualen+1));
      ua2a_map[0] = -1;
      for(apos = 1; apos <= msa->alen; apos++)
	if(a2ua_map[apos] != -1)
	  ua2a_map[a2ua_map[apos]] = apos;

      ESL_ALLOC(uasubseq, sizeof(char) * (sq[i]->n+1));
      esl_abc_Textize(msa->abc, sq[i]->dsq, sq[i]->n, uasubseq);
      esl_abc_Textize(msa->abc, msa->ax[i], msa->alen, aseq);

      esl_strdup(aseq, -1, &(uaseq));
      esl_strdealign(uaseq, uaseq, "-_.", NULL);
      offset = strstr(uaseq, uasubseq);
      if(offset == NULL) ESL_XFAIL(eslEINVAL, errbuf, "in trim_msa(), sq[%d] is not a subseq of msa seq %d\n", i, i);
      uastart = offset  - uaseq + 1;
      uaend   = uastart + strlen(uasubseq) - 1;
      astart  = ua2a_map[uastart];
      aend    = ua2a_map[uaend];
      free(ua2a_map);
      free(a2ua_map);

      for(apos = 1;        apos <  astart;    apos++) msa->ax[i][apos] = msa->abc->K; /* make it a gap */
      for(apos = aend + 1; apos <= msa->alen; apos++) msa->ax[i][apos] = msa->abc->K; /* make it a gap */
      free(uaseq);
      free(uasubseq);
    }

  for(i = 0; i < msa->nseq; i++)
    esl_sq_Destroy(sq[i]);
  free(sq);
      
  free(aseq);
  return eslOK;

 ERROR:
  return status;
}

/* dump_insert_info
 *                   
 * Given an MSA with RF annotation, print out information about how many 'insertions' come
 * after each non-gap RF column (consensus column). 
 */
static int dump_insert_info(FILE *fp, ESL_MSA *msa, char *errbuf)
{
  int status;
  int apos, cpos;
  int **ict;
  int *total_ict, *med_ict;
  int i, l;
  int clen;
  int nseq;
  int *len;

  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL)) ESL_XFAIL(eslEINVAL, errbuf, "in dump_insert_info(), msa must be digitized.");
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --iplot.");

  ESL_ALLOC(total_ict,  sizeof(int) * (msa->alen+2));
  ESL_ALLOC(med_ict,  sizeof(int) * (msa->alen+2));
  esl_vec_ISet(total_ict, (msa->alen+2), 0);
  esl_vec_ISet(med_ict, (msa->alen+2), 0);

  ESL_ALLOC(ict,  sizeof(int *) * (msa->alen+2));
  for(i = 0; i <= msa->alen; i++)
    {
      ESL_ALLOC(ict[i],  sizeof(int) * (msa->nseq));
      esl_vec_ISet(ict[i], (msa->nseq), 0);
    }

  fprintf(fp, "# %8s  %10s  %8s  %8s  %8s\n", "cons col", "nseq w/ins",  "freq ins", "avg len",  "med len");
  fprintf(fp, "# %8s  %10s  %8s  %8s  %8s\n", "--------", "----------", "--------", "--------", "--------");

  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++)
    {
      if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) 
	cpos++;
      else
	for(i = 0; i < msa->nseq; i++)
	  if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) { 
	    ict[cpos][i]++;
	    total_ict[cpos]++;
	  }	  
    }
  clen = cpos;

  /* determine avg median length for each insertion */
  for(cpos = 0; cpos <= clen; cpos++)
    {
      if(total_ict[cpos] > 0) { 
	nseq = 0;
	for(i = 0; i < msa->nseq; i++) { 
	  if(ict[cpos][i] >= 1) nseq++;
	}
	ESL_ALLOC(len, sizeof(int) * nseq);
	l = 0;
	for(i = 0; i < msa->nseq; i++) { 
	  if(ict[cpos][i] >= 1)
	    len[l++] = ict[cpos][i];
	}
	qsort(len, nseq, sizeof(int), compare_ints);
	med_ict[cpos] = len[nseq / 2];
	free(len);
      }      
    }
  for(cpos = 0; cpos <= clen; cpos++)
    {
      nseq = 0;
      for(i = 0; i < msa->nseq; i++) if(ict[cpos][i] >= 1) nseq++;
      if(nseq > 0) 
	fprintf(fp, "  %8d  %10d  %8.6f  %8.3f  %8d\n", cpos, nseq, (float) nseq / (float) msa->nseq, ((float) total_ict[cpos] / (float) nseq), med_ict[cpos]);
    }

  for(i = 0; i <= msa->alen; i++)
    free(ict[i]);
  free(ict);
  free(total_ict);
  free(med_ict);

  return eslOK;

 ERROR:
  return status;
}


/* dump_residue_info
 *                   
 * Given an MSA, print out the number of sequences with
 * a non-gap residue in each column of the alignment.
 */
static int dump_residue_info(FILE *fp, ESL_MSA *msa, char *errbuf)
{
  int status;
  int apos, cpos;
  int rct;
  int i;
  int has_rf;

  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL)) ESL_XFAIL(eslEINVAL, errbuf, "in dump_residue_info(), msa must be digitized.");
  has_rf = (msa->rf == NULL) ? FALSE : TRUE;

  if(has_rf) { 
    fprintf(fp, "# %8s  %7s  %8s  %8s\n", "cons col", "aln col", "num res",  "freq res");
    fprintf(fp, "# %8s  %7s  %8s  %8s\n", "--------", "-------", "--------", "--------");
  }  
  else { 
    fprintf(fp, "# %7s  %8s  %8s\n", "aln col", "num res",  "freq res");
    fprintf(fp, "# %7s  %8s  %8s\n", "-------", "--------", "--------");
  }
  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++) {
    rct = 0;
    if(has_rf && (! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)]))) cpos++;
    for(i = 0; i < msa->nseq; i++) { 
      if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) rct++; 
    }

    if(has_rf) fprintf(fp, "  %8d  %7d  %8d  %8.6f\n", cpos, apos, rct, (float) rct / (float) msa->nseq);
    else       fprintf(fp, "  %7d  %8d  %8.6f\n", apos, rct, (float) rct / (float) msa->nseq);
  }

  return eslOK;

 ERROR:
  return status;
}


/* dump_delete_info
 *                   
 * Given an MSA, print out the number of sequences with
 * gaps residue in each consensus column of the alignment.
 */
static int dump_delete_info(FILE *fp, ESL_MSA *msa, char *errbuf)
{
  int status;
  int apos, cpos;
  int dct;
  int i;

  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL)) ESL_XFAIL(eslEINVAL, errbuf, "in dump_residue_info(), msa must be digitized.");
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --dinfo.");

  fprintf(fp, "# Number of sequences in file: %d\n", msa->nseq);
  fprintf(fp, "# Only non-gap RF columns with > 0 deletes are listed.\n");
  fprintf(fp, "#\n");
  fprintf(fp, "# %8s  %7s  %8s  %8s\n", "cons col", "aln col", "num del",  "freq del");
  fprintf(fp, "# %8s  %7s  %8s  %8s\n", "--------", "-------", "--------", "--------");
  
  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++) {
    dct = 0;
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) { 
      cpos++;
      for(i = 0; i < msa->nseq; i++) { 
	if(esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) dct++; 
      }
      if(dct > 0) fprintf(fp, "  %8d  %7d  %8d  %8.6f\n", cpos, apos, dct, (float) dct / (float) msa->nseq);
    }
  }
  return eslOK;

 ERROR:
  return status;
}

/* plot_inserts
 *                   
 * Given an MSA with RF annotation, print a postscript heatmap of how
 * many insertions are after each non-gap RF column (consensus column)
 * in each sequence.
 */
static int plot_inserts(FILE *fp, ESL_MSA *msa, int do_log, char *errbuf)
{
  int status;
  int apos, cpos;
  int i;
  int clen;
  ESL_DMATRIX *I;

  /* contract check */
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --iplot.");
  if(! (msa->flags & eslMSA_DIGITAL))
    ESL_XFAIL(eslEINVAL, errbuf, "in plot_inserts(), msa must be digitized.");

  clen = 0;
  for(apos = 1; apos <= msa->alen; apos++)
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) clen++;

  I = esl_dmatrix_Create(msa->nseq, (clen+1));
  esl_dmatrix_SetZero(I);

  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++) {
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) 
      cpos++;
    else
      for(i = 0; i < msa->nseq; i++)
	if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) I->mx[i][cpos] += 1.;
  }
  if(do_log) {
    for(i = 0; i < msa->nseq; i++)
      for(cpos = 0; cpos <= clen; cpos++)
	if(I->mx[i][cpos] > 0) 
	  I->mx[i][cpos] = log(I->mx[i][cpos]);
	else 
	  I->mx[i][cpos] = -1; /* don't want 0s to change to -inf, for coloring scheme */
  }
  else { 
    for(i = 0; i < msa->nseq; i++)
      for(cpos = 0; cpos <= clen; cpos++)
	if(I->mx[i][cpos] == 0) I->mx[i][cpos] = (-1 * esl_dmx_Max(I)) / 2; /* for better resolution on heatmap */
  }

  /* dmx_Visualize(fp, I, esl_dmx_Min(I), esl_dmx_Max(I)); */
  dmx_Visualize(fp, I, (-1 * esl_dmx_Max(I)), esl_dmx_Max(I));
  /* esl_dmatrix_Dump(stdout, I, NULL, NULL); */
  esl_dmatrix_Destroy(I);
  return eslOK;

 ERROR:
  return status;
}


/* plot_gaps
 *                   
 * Given an MSA with RF annotation, print a postscript checkboard grid 
 * showing which sequences have gaps in each non-gap RF column. 
 */
static int plot_gaps(FILE *fp, ESL_MSA *msa, char *errbuf)
{
  int status;
  int apos, cpos;
  int i;
  int clen;
  ESL_DMATRIX *G;

  /* contract check */
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --gplot.");
  if(! (msa->flags & eslMSA_DIGITAL))
    ESL_XFAIL(eslEINVAL, errbuf, "in plot_gaps(), msa must be digitized.");

  clen = 0;
  for(apos = 1; apos <= msa->alen; apos++)
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) clen++;

  G = esl_dmatrix_Create(msa->nseq, (clen+1));
  esl_dmatrix_SetZero(G);

  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++) {
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) {
      cpos++;
      for(i = 0; i < msa->nseq; i++)
	if(esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) G->mx[i][cpos] += 1.;
    }
  }
  /* dmx_Visualize(fp, G, esl_dmx_Min(G), esl_dmx_Max(G)); */
  dmx_Visualize(fp, G, -1, 1);
  esl_dmatrix_Destroy(G);
  return eslOK;

 ERROR:
  return status;
}

/* get_tree_order
 *                   
 * Given a tree, determine the branching order of the sequences
 * it represents by traversing it preorder.
 */
static int get_tree_order(ESL_TREE *T, char *errbuf, int **ret_order)
{
  int status;
  int opos = 0;
  int nd;
  int *order; 
  ESL_STACK *pda;
  ESL_ALLOC(order, sizeof(int) * T->N);

  opos = 0;
  pda  = esl_stack_ICreate();
  esl_stack_IPush(pda, T->right[0]);
  esl_stack_IPush(pda, T->left[0]);
  while (esl_stack_IPop(pda, &nd) != eslEOD)
    {
      if (nd > 0) { /* a node */
	esl_stack_IPush(pda, T->right[nd]); /* index for right child */
	esl_stack_IPush(pda, T->left[nd]);  /* index for left child */
      }
      else /* nd <= 0, a child */
	order[opos++] = nd * -1;
    }
  *ret_order = order;
  esl_stack_Destroy(pda);
  return eslOK;

 ERROR:
  return status;
}


/* reorder_msa
 *                   
 * Given an array specifying a new order for the sequences in
 * the MSA, reorder it by swapping pointers.
 */
static int
reorder_msa(ESL_MSA *msa, int *order, char *errbuf)
{
  int status;
  char **tmp; 
  ESL_ALLOC(tmp, sizeof(char *) * msa->nseq);
  int i, a;

  /* contract check */
  /* 'order' must be have nseq elements, elements must be in range [0..nseq-1], no duplicates  */
  int *covered;
  ESL_ALLOC(covered, sizeof(int) * msa->nseq);
  esl_vec_ISet(covered, msa->nseq, 0);
  for(i = 0; i < msa->nseq; i++) { 
    if(covered[order[i]]) ESL_FAIL(eslEINVAL, errbuf, "reorder_msa() order array has duplicate entries for i: %d\n", i);
    covered[order[i]] = 1;
  }
  free(covered);

  /* swap aseq or ax (one or the other must be non-NULL) */
  if(msa->flags & eslMSA_DIGITAL) { /* digital MSA */
    ESL_DSQ **tmp_dsq; 
    ESL_ALLOC(tmp_dsq, sizeof(ESL_DSQ *) * msa->nseq);
    for(i = 0; i < msa->nseq; i++) tmp_dsq[i] = msa->ax[i];
    for(i = 0; i < msa->nseq; i++) msa->ax[i] = tmp_dsq[order[i]];
    free(tmp_dsq);
  }
  else { /* text MSA */
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->aseq[i];
    for(i = 0; i < msa->nseq; i++) msa->aseq[i] = tmp[order[i]];
  }

  /* swap sqnames (mandatory) */
  for(i = 0; i < msa->nseq; i++) tmp[i] = msa->sqname[i];
  for(i = 0; i < msa->nseq; i++) msa->sqname[i] = tmp[order[i]];

  /* swap sqacc, if they exist */
  if(msa->sqacc != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->sqacc[i];
    for(i = 0; i < msa->nseq; i++) msa->sqacc[i] = tmp[order[i]];
  }

  /* swap sqdesc, if they exist */
  if(msa->sqdesc != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->sqdesc[i];
    for(i = 0; i < msa->nseq; i++) msa->sqdesc[i] = tmp[order[i]];
  }

  /* swap ss, if they exist */
  if(msa->ss != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->ss[i];
    for(i = 0; i < msa->nseq; i++) msa->ss[i] = tmp[order[i]];
  }

  /* swap sa, if they exist */
  if(msa->sa != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->sa[i];
    for(i = 0; i < msa->nseq; i++) msa->sa[i] = tmp[order[i]];
  }

  /* swap gs annotation, if it exists */
  for(a = 0; a < msa->ngs; a++) {
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->gs[a][i];
    for(i = 0; i < msa->nseq; i++) msa->gs[a][i] = tmp[order[i]];
  }

  /* swap gr annotation, if it exists */
  for(a = 0; a < msa->ngr; a++) {
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->gr[a][i];
    for(i = 0; i < msa->nseq; i++) msa->gr[a][i] = tmp[order[i]];
  }
  free(tmp);
  return eslOK;

 ERROR: 
  return status;
}

/****************************************************************
 * Stolen from hmmer/h3/heatmap.c SVN revision 2171
 * as dmx_Visualize. Then modified so that the full 
 * matrix is printed (not half split diagonally).
 */
/* dmx_Visualize()
 * Incept:    SRE, Wed Jan 24 11:58:21 2007 [Janelia]
 *
 * Purpose:   
 *            
 *            Color scheme roughly follows Tufte, Envisioning
 *            Information, p.91, where he shows a beautiful
 *            bathymetric chart. The CMYK values conjoin two
 *            recommendations from ColorBrewer (Cindy Brewer
 *            and Mark Harrower) 
 *            [http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer.html],
 *            specifically the 9-class sequential2 Blues and
 *            9-class sequential YlOrBr.
 * 
 *            Might eventually become part of Easel, once mature?
 *           
 * Note:      Binning rules basically follow same convention as
 *            esl_histogram. nb = xmax-xmin/w, so w = xmax-xmin/nb; 
 *            picking bin is (int) ceil((x - xmin)/w) - 1. (xref
 *            esl_histogram_Score2Bin()). This makes bin b contain
 *            values bw+min < x <= (b+1)w+min. (Which means that 
 *            min itself falls in bin -1, whoops - but we catch
 *            all bin<0 and bin>=nshades and put them in the extremes.
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
dmx_Visualize(FILE *fp, ESL_DMATRIX *D, double min, double max)
 {
   int    nshades   = 18;
   double cyan[]    = { 1.00, 1.00, 0.90, 0.75, 0.57, 0.38, 0.24, 0.13, 0.03,
			0.00, 0.00, 0.00, 0.00, 0.00, 0.07, 0.20, 0.40, 0.60};
   double magenta[] = { 0.55, 0.45, 0.34, 0.22, 0.14, 0.08, 0.06, 0.03, 0.01,
			0.00, 0.03, 0.11, 0.23, 0.40, 0.55, 0.67, 0.75, 0.80};
   double yellow[]  = { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
			0.10, 0.25, 0.40, 0.65, 0.80, 0.90, 1.00, 1.00, 1.00};
   double black[]   = { 0.30, 0.07, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
			0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
   double w;			
   int    i,j;
   int    bin;
   int    boxsize;		/* box size in points */
   int    xcoord, ycoord;	/* postscript coords in points */
   int    leftmargin, rightmargin;
   int    bottommargin, topmargin;
   float  fboxsize;		/* box size in fractional points */

   /* Set some defaults that might become arguments later.
    */
   leftmargin   = rightmargin = 20;
   bottommargin = topmargin   = 20;

   /* Determine some working parameters 
    */
   w = (max-min) / (double) nshades; /* w = bin size for assigning values->colors*/
   boxsize = ESL_MAX(1, (ESL_MIN((792 - bottommargin) / D->n, 
				 (612 - leftmargin)   / D->m)));
   fboxsize= ESL_MIN( (792. - ((float) bottommargin + topmargin))   / (float) D->n, 
		      (612. - ((float) leftmargin   + rightmargin)) / (float) D->m);


   fprintf(fp, "%.4f %.4f scale\n", (fboxsize/(float) boxsize), (fboxsize/(float) boxsize));
   /* printf("n: %d\nm: %d\n", D->n, D->m); */
   for (i = 0; i < D->n; i++) {
     /* printf("\n"); */
     /* for (j = i; j < D->n; j++) */
     for (j = 0; j < D->m; j++)
       {
	 /* printf("i: %4d j: %4d %5.1f\n", i, j, D->mx[i][j]); */
	 xcoord = j * boxsize + leftmargin;
	 ycoord = (D->m-(i+1)) * boxsize + bottommargin; /* difference w/heatmap.c: (D->m-i+1) */
	 
	 if      (D->mx[i][j] == -eslINFINITY) bin = 0;
	 else if (D->mx[i][j] ==  eslINFINITY) bin = nshades-1;
	 else {
	   bin    = (int) ceil((D->mx[i][j] - min) / w) - 1;
	   if (bin < 0)        bin = 0;
	   if (bin >= nshades) bin = nshades-1;
	 }
	 
	fprintf(fp, "newpath\n");
	fprintf(fp, "  %d %d moveto\n", xcoord, ycoord);
	fprintf(fp, "  0  %d rlineto\n", boxsize);
	fprintf(fp, "  %d 0  rlineto\n", boxsize);
	fprintf(fp, "  0 -%d rlineto\n", boxsize);
	fprintf(fp, "  closepath\n");
  	fprintf(fp, " %.2f %.2f %.2f %.2f setcmykcolor\n",
		cyan[bin], magenta[bin], yellow[bin], black[bin]);
	fprintf(fp, "  fill\n");
      }
   }
  fprintf(fp, "showpage\n");
  return eslOK;
}


/* read_mask_file
 *
 * Given an open file pointer, read the first token of the
 * file and return it as *ret_mask. It must contain only
 * '0' or '1' characters.
 *
 * Returns:  eslOK on success.
 */
int
read_mask_file(char *filename, char *errbuf, char **ret_mask, int *ret_mask_len)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  char           *mask;
  int             toklen;
  int             n;

  if (esl_fileparser_Open(filename, &efp) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to open %s in read_mask_file\n", filename);
  esl_fileparser_SetCommentChar(efp, '#');
  
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to read a single token from %s\n", filename);

  ESL_ALLOC(mask, sizeof(char) * (toklen+1));
  for(n = 0; n < toklen; n++) { 
    if((tok[n] == '0') || (tok[n] == '1')) { 
      mask[n] = tok[n];
    }
    else { ESL_FAIL(eslFAIL, errbuf, "read a non-0 and non-1 character (%c) in the mask file %s\n", tok[n], filename); }
  }
  mask[n] = '\0';

  *ret_mask = mask;
  *ret_mask_len = n;
  esl_fileparser_Close(efp);
  return eslOK;
  
 ERROR:
  return eslEMEM;
}


/* handle_post_opts
 *                   
 * Read "#=GR POST" annotation into a 2D matrix, each sequence
 * is a row, each residue is a column. Handle any command line
 * options that use the posterior info.
 *
 */      
static int handle_post_opts(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa)
{
  int    status;
  int    s,c;                 /* counters over sequences, columns of MSA */
  int   *nongap_c, *nongap_s; /* number of non-gap posterior values for each column/sequence respectively */
  float *sum_c, *sum_s;       /* sum of non-gap posterior values for each column/sequence respectively */
  float *min_c, *min_s;       /* min of non-gap posterior values for each column/sequence respectively */
  float *avg_c, *avg_s;       /* average non-gap posterior values for each column/sequence respectively */
  int   *athresh_c;           /* [0..c..msa->alen-1] number of sequences with residue in this column with post value >= pthresh */
  float *athresh_fract_c;     /* [0..c..msa->alen-1] fraction of non-gap sequences with residue in this column with post value >= pthresh */
  int    ridx1, ridx2;
  int    r;
  float  p;
  int    do_pfract = (! esl_opt_IsDefault(go, "--pfract"));
  int    do_prf    = (! esl_opt_IsDefault(go, "--p-rf"));
  int    do_pinfo  = (! esl_opt_IsDefault(go, "--pinfo"));
  float  pfract;
  float  pthresh   =    esl_opt_GetReal(go, "--pthresh"); /* default is 0.95 */
  int  *useme; 
  int *c2a_map;
  int clen;
  int cpos;
  int nkept;
  int ndigits;
  int ir1, ir2;
  int nongap_total = 0;
  int nongap_total_rf = 0;
  float sum_total = 0.;
  float sum_total_rf = 0.;
  FILE *pinfofp = NULL;  /* output file for --pinfo */

  if((!do_pfract) && (!do_pinfo)) ESL_FAIL(eslEINVAL, errbuf, "handle_post_opts(): --pinfo nor --pfract options selected, shouldn't be in this function.");

  /* Find out which #=GR line is the POST, Post, or post line (if more than one exist, last one is chosen) */
  ridx1 = -1;
  ridx2 = -1;
  ndigits = 0;
  for (r = 0; r < msa->ngr; r++) { 
    if (strcmp(msa->gr_tag[r], "POST")   == 0) { ridx1 = r; ndigits = 1; }
    if (strcmp(msa->gr_tag[r], "Post")   == 0) { ridx1 = r; ndigits = 1; }
    if (strcmp(msa->gr_tag[r], "post")   == 0) { ridx1 = r; ndigits = 1; }
    if (strcmp(msa->gr_tag[r], "POSTX.") == 0) { ridx1 = r; ndigits = 1; }
    if (strcmp(msa->gr_tag[r], "POST.X") == 0) { ridx2 = r; ndigits = 2; }
  }
  if(ndigits == 1 && ridx1 == -1) { 
    if(do_pfract) ESL_FAIL(eslEINVAL, errbuf, "--pfract requires \"#=GR POST\", \"#=GR Post\", \"#=GR post\", \"#=GR POSTX.\", or \"#=GR POSTX.\" and \"#=GR POST.X\" annotation in %s.\n", esl_opt_GetArg(go,1));
    if(do_pinfo)  ESL_FAIL(eslEINVAL, errbuf, "--pinfo  requires \"#=GR POST\", \"#=GR Post\", \"#=GR post\", \"#=GR POSTX.\", or \"#=GR POSTX.\" and \"#=GR POST.X\" annotation in %s.\n", esl_opt_GetArg(go,1));
  }
  if(ndigits == 2 && (ridx1 == -1 || ridx2 == -1)) { 
    if(do_pfract) ESL_FAIL(eslEINVAL, errbuf, "--pfract requires \"#=GR POST\", \"#=GR Post\", \"#=GR post\", or \"#=GR POSTX.\" and \"#=GR POST.X\" annotation in %s.\n", esl_opt_GetArg(go,1));
    if(do_pinfo)  ESL_FAIL(eslEINVAL, errbuf, "--pinfo  requires \"#=GR POST\", \"#=GR Post\", \"#=GR post\", or \"#=GR POSTX.\" and \"#=GR POST.X\" annotation in %s.\n", esl_opt_GetArg(go,1));
  }
  if(msa->rf == NULL) { 
    if(do_prf)    ESL_FAIL(eslEINVAL, errbuf, "--p-rf requires \"#=GC RF\" annotation in %s.\n", esl_opt_GetArg(go,1));
  }
  /*ESL_ALLOC(post_sc, sizeof(int *) * msa->nseq);
    for(i = 0; i < msa->nseq; i++) 
    ESL_ALLOC(post_sc, sizeof(int) * msa->alen);
    ESL_ALLOC(post_cs, sizeof(int *) * msa->alen);
    for(i = 0; i < msa->alen; i++) 
    ESL_ALLOC(post_cs, sizeof(int) * msa->nseq);
  */
     
  ESL_ALLOC(nongap_c, sizeof(int) * msa->alen);
  ESL_ALLOC(sum_c,    sizeof(float) * msa->alen);
  ESL_ALLOC(min_c,    sizeof(float) * msa->alen);
  ESL_ALLOC(avg_c,    sizeof(float) * msa->alen);
  ESL_ALLOC(athresh_c,sizeof(int)   * msa->alen);
  ESL_ALLOC(athresh_fract_c,sizeof(float) * msa->alen);
  esl_vec_ISet(nongap_c, msa->alen, 0);
  esl_vec_FSet(sum_c,    msa->alen, 0.);
  esl_vec_FSet(min_c,    msa->alen, 10.);
  esl_vec_ISet(athresh_c,msa->alen, 0);

  ESL_ALLOC(nongap_s, sizeof(int) * msa->nseq);
  ESL_ALLOC(sum_s, sizeof(float) * msa->nseq);
  ESL_ALLOC(min_s, sizeof(float) * msa->nseq);
  ESL_ALLOC(avg_s, sizeof(float) * msa->nseq);
  esl_vec_ISet(nongap_s, msa->nseq, 0);
  esl_vec_FSet(sum_s,    msa->nseq, 0.);
  esl_vec_FSet(min_s,    msa->nseq, 10.);

  if(ndigits == 1) {    
    for(s = 0; s < msa->nseq; s++) { 
      for(c = 0; c < msa->alen; c++) { 
	if(! esl_abc_CIsGap(msa->abc, msa->gr[ridx1][s][c])) {
	  switch(msa->gr[ridx1][s][c]) { 
	  case '*': p = 1.0; break;
	  case '9': p = 0.9; break;
	  case '8': p = 0.8; break;
	  case '7': p = 0.7; break;
	  case '6': p = 0.6; break;
	  case '5': p = 0.5; break;
	  case '4': p = 0.4; break;
	  case '3': p = 0.3; break;
	  case '2': p = 0.2; break;
	  case '1': p = 0.1; break;
	  case '0': p = 0.0; break;
	  default: 
	    ESL_FAIL(eslEINVAL, errbuf, "reading post annotation for seq: %d aln column: %d, unrecognized residue: %c\n", s, c, msa->gr[ridx1][s][c]);
	  }
	  sum_c[c] += p;
	  sum_s[s] += p;
	  nongap_c[c]++;
	  nongap_s[s]++;
	  min_c[c] = ESL_MIN(min_c[c], p);
	  min_s[s] = ESL_MIN(min_s[s], p);
	  if(p >= pthresh) athresh_c[c]++;
	}
	else p = -1; /* gap */
	
	/*post_sc[s][c] = p;
	  post_cs[c][s] = p;*/
      }
    }
  }
  if(ndigits == 2) { 
    for(s = 0; s < msa->nseq; s++) { 
      for(c = 0; c < msa->alen; c++) { 
	if(! esl_abc_CIsGap(msa->abc, msa->gr[ridx1][s][c])) {
	  if(esl_abc_CIsGap(msa->abc, msa->gr[ridx2][s][c])) ESL_FAIL(eslEINVAL, errbuf, "reading post annotation for seq: %d aln column: %d, post 'tens' value non-gap but post 'ones' value is gap.\n", s, c);
	  if(msa->gr[ridx1][s][c] == '*') {
	    if(msa->gr[ridx2][s][c] != '*') ESL_FAIL(eslEINVAL, errbuf, "reading post annotation for seq: %d aln column: %d, post 'tens' value '*' but post 'ones' value != '*'.\n", s, c);
	    p = 1.0;
	  }
	  else {
	    ir1 = (int) (msa->gr[ridx1][s][c] - '0');
	    ir2 = (int) (msa->gr[ridx2][s][c] - '0');
	    p = ((float) ir1 * 10. + ir2) * .01;
	    /* printf("r1: %c %d r2: %c %d p: %.2f\n", msa->gr[ridx1][s][c], ir1, msa->gr[ridx2][s][c], ir2, p);*/
	  }
	  sum_c[c] += p;
	  sum_s[s] += p;
	  nongap_c[c]++;
	  nongap_s[s]++;
	  min_c[c] = ESL_MIN(min_c[c], p);
	  min_s[s] = ESL_MIN(min_s[s], p);
	  if(p >= pthresh) athresh_c[c]++;
	}
	else p = -1; /* gap */
	
	/*post_sc[s][c] = p;
	  post_cs[c][s] = p;*/
      }
    }
  }

  if(msa->rf != NULL) map_cpos_to_apos(msa, &c2a_map, &clen);
  else c2a_map = NULL;

  /* get averages */
  for(s = 0; s < msa->nseq; s++) { 
    avg_s[s]  =  (float) sum_s[s] / (float) nongap_s[s];
  }
  cpos = 1;
  for(c = 0; c < msa->alen; c++) { 
    avg_c[c]  = (float) sum_c[c] / (float) nongap_c[c];
    sum_total += sum_c[c];
    nongap_total += nongap_c[c];
    if(c2a_map != NULL) {
      if(c2a_map[cpos] == (c+1)) {  /* off-by-one, c2a_map is 1..clen, c is 0..alen */
	cpos++; 
	sum_total_rf += sum_c[c];
	nongap_total_rf += nongap_c[c];
      }
    }
  }
  /* determine the fraction of sequences in each column with POSTs that exceed pthresh */
  for(c = 0; c < msa->alen; c++) 
    athresh_fract_c[c] = (nongap_c[c] > 0) ? ((float) athresh_c[c] / (float) nongap_c[c]) : 0;


  printf("\nAverage posterior value:                            %.5f (%d non-gap residues)\n", (float) sum_total / (float) nongap_total, nongap_total);
  if(c2a_map != NULL) 
    printf("Average posterior value in non-gap #=GC RF columns: %.5f (%d non-gap RF residues)\n", (float) sum_total_rf / (float) nongap_total_rf, nongap_total_rf);
  printf("\n");


  /* if nec, print posterior info */
  cpos = 1;
  if(do_pinfo) { 
    if ((pinfofp = fopen(esl_opt_GetString(go, "--pinfo"), "w")) == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open --pinfo output file %s\n", esl_opt_GetString(go, "--pinfo"));
    fprintf(pinfofp, "# Posterior stats per column:\n");
    if(msa->rf != NULL) { 
      fprintf(pinfofp, "# %3s %5s %6s %6s %6s > %5.3f\n", "rf?", "col", "nongap", "avg", "min", pthresh);
      fprintf(pinfofp, "# %3s %5s %6s %6s %6s %7s\n", "---", "-----", "------", "------", "------", "-------");
      for(c = 0; c < msa->alen; c++) { 
	if(c2a_map[cpos] == (c+1)) { /* off-by-one, c2a_map is 1..clen, c is 0..alen */
	  cpos++; 
	  fprintf(pinfofp, "  *   "); 
	}  
	else fprintf(pinfofp, "      ");
	if(nongap_c[c] == 0) fprintf(pinfofp, "%5d %6.3f %6.3f %6.1f %7.3f\n", c+1, ((float) (nongap_c[c]) / ((float) msa->nseq)),      0.0,      0.0, athresh_fract_c[c]);
	else                 fprintf(pinfofp, "%5d %6.3f %6.3f %6.1f %7.3f\n", c+1, ((float) (nongap_c[c]) / ((float) msa->nseq)), avg_c[c], min_c[c], athresh_fract_c[c]);
      }
    }
    else { /* msa->rf is NULL, we can't indicate the non-gap RF columns */
      fprintf(pinfofp, "%5s %6s %6s %6s > %5.3f\n", "col", "nongap", "avg", "min", pthresh);
      fprintf(pinfofp, "%5s %6s %6s %6s %7s\n", "-----", "------", "------", "------", "-------");
      for(c = 0; c < msa->alen; c++) 
	fprintf(pinfofp, "%5d %6.3f %6.3f %6.1f %7.3f\n", c+1, ((float) (nongap_c[c]) / ((float) msa->nseq)), avg_c[c], min_c[c], athresh_fract_c[c]);
    }
    fprintf(pinfofp, "\n\n");

    fprintf(pinfofp, "# Posterior stats per sequence:\n");
    fprintf(pinfofp, "# %5s %-60s %6s %6s %6s\n", "idx",   "seq name", "nongap", "avg", "min");
    fprintf(pinfofp, "# %5s %-60s %6s %6s %6s\n", "-----", "------------------------------------------------------------", "------", "------", "------");
    for(s = 0; s < msa->nseq; s++) { 
      fprintf(pinfofp, "  %5d %-60s %6.3f %6.3f %6.2f\n", s+1, msa->sqname[s], ((float) (nongap_s[s]) / ((float) msa->alen)), avg_s[s], min_s[s]); 
    }
    fclose(pinfofp);
  }

  /* optionally, add/rewrite msa->rf if --pfract enabled */
  if(do_pfract) { 
    pfract = esl_opt_GetReal(go, "--pfract");  
    ESL_ALLOC(useme, sizeof(int) * (msa->alen+1));
    if(do_prf) { /* only look at consensus columns */
      cpos = 1;
      for(c = 0; c < msa->alen; c++) { 
	if(c2a_map[cpos] == (c+1)) { /* off-by-one, c2a_map is 1..clen, c is 0..alen */
	  cpos++;
	  if(athresh_fract_c[c] >= pfract) useme[c] = 1; 
	  else                             useme[c] = 0; 
	}
	else useme[c] = 0; 
      }    
    }
    else { /* look at all columns */
      for(c = 0; c < msa->alen; c++) {  
	if(athresh_fract_c[c] >= pfract) useme[c] = 1; 
	else                             useme[c] = 0; 
      }
    }
    useme[msa->alen] = '\0';

    nkept = 0;
    write_rf_given_useme(go, errbuf, msa, useme);
    for(c = 0; c < msa->alen; c++) nkept += useme[c];
    free(useme);
    if(do_prf)  printf("\n%d of %d RF columns (%.3f) pass threshold\n\n", nkept, clen, (float) nkept / (float) clen);
    else        printf("\n%d of %" PRId64 " columns (%.3f) pass threshold\n\n", nkept, msa->alen, (float) nkept / (float) msa->alen);
  }

  /*for(s = 0; s < msa->nseq; s++) free(post_sc[s]);
    free(post_sc);
    for(c = 0; c < msa->alen; c++) free(post_cs[c]);
    free(post_cs);*/
  
  free(athresh_fract_c);
  free(athresh_c);
  free(nongap_s);
  free(nongap_c);
  free(min_c);
  free(min_s);
  free(avg_c);
  free(avg_s);
  free(sum_s);
  free(sum_c);
  if(c2a_map != NULL) free(c2a_map);

  return eslOK;

 ERROR:
  return status;
}


/* output_rf_as_mask
 *
 * Given an MSA with rf annotation, convert it to a lanemask of 1s and 0s.
 * 1s for non-gap RF columns, 0s for gap RF columns.
 */
static int
output_rf_as_mask(FILE *fp, char *errbuf, ESL_MSA *msa)
{
  int status;
  int  apos;
  char *mask;

  if(msa->rf == NULL) ESL_FAIL(eslEINVAL, errbuf, "msa->rf is NULL, and we're trying to convert it to a 1/0 mask.");
  ESL_ALLOC(mask, sizeof(char) * (msa->alen+1));

  for (apos = 0; apos < msa->alen; apos++) 
    if(esl_abc_CIsGap(msa->abc, msa->rf[apos]))  mask[apos] = '0';
    else                                         mask[apos] = '1';
  mask[msa->alen] = '\0';

  fprintf(fp, "%s\n", mask);
  free(mask);

  return eslOK;
 ERROR:
  return status;
}

/* expand_msa2mask
 *
 * Given an MSA <msa> and a lanemask <xmask> with exactly msa->alen 1s in it.
 * Add 100% gap columns in between each column as dictated by <xmask>.
 *
 * For example if lanemask is 100101, msa->alen is 3, we add 2 100% gap
 * columns after column 1, and 1 100% gap column after column 2, to make
 * the msa length = length(xmask) = 6.
 */
static int
expand_msa2mask(char *errbuf, ESL_MSA *msa, char *xmask, ESL_MSA **newmsa)
{
  int status;
  int  mpos;
  int  masklen;
  int *nzeroesA;
  int  nones = 0;

  if(xmask == NULL) ESL_FAIL(eslEINVAL, errbuf, "expand_msa2mask(), xmask is NULL.");

  masklen = strlen(xmask);
  /* count 1s in xmask */
  for (mpos = 0; mpos < masklen; mpos++) { 
    if     (xmask[mpos] == '1') nones++;
    else if(xmask[mpos] == '0') ; /* do nothing */
    else    ESL_FAIL(eslEINVAL, errbuf, "--xmask mask char number %d is not a 1 nor a 0, but a %c\n", mpos+1, xmask[mpos]);
  }
  if(nones != msa->alen) ESL_FAIL(eslEINVAL, errbuf, "expand_msa2mask(), number of 1s in --xmask file: %d != msa->alen: %" PRId64 ", they must be equal.", nones, msa->alen);

  /* determine number of 0s after each consensus column */
  nones = 0;
  ESL_ALLOC(nzeroesA, sizeof(int) * masklen+1);
  esl_vec_ISet(nzeroesA, (masklen+1), 0);
  for (mpos = 0; mpos < masklen; mpos++) { 
    if     (xmask[mpos] == '1') nones++;
    else if(xmask[mpos] == '0') nzeroesA[nones]++;
    else    ESL_FAIL(eslEINVAL, errbuf, "--xmask mask char number %d is not a 1 nor a 0, but a %c\n", mpos+1, xmask[mpos]);
  }
  
  /*int i;
  for (i = 0; i <= nones; i++) { 
    printf("nzeroes[%3d]: %3d\n", i, nzeroesA[i]);
    }*/

  /* add the 100% gap columns */
  if((status = add_gap_columns_to_msa(errbuf, msa, nzeroesA, newmsa, TRUE)) != eslOK) return status ;
  /* new alen should equal masklen */
  if((*newmsa)->alen != masklen) ESL_FAIL(eslEINVAL, errbuf, "expand_msa2mask(), new msa->alen: (%" PRId64 ") != length of mask (%d), this shouldn't happen.", (*newmsa)->alen, masklen);
  free(nzeroesA);

  return eslOK;
 ERROR:
  return status;
}

/* Function: compare_ints()
 * 
 * Purpose:  Comparison function for qsort(). Used 
 *           by msa_median_length().
 */ 
static int 
compare_ints(const void *el1, const void *el2)
{
  if      ((* ((int *) el1)) > (* ((int *) el2)))  return 1;
  else if ((* ((int *) el1)) < (* ((int *) el2)))  return 1;
  return 0;
}

/* Function: msa_median_length()
 * 
 * Purpose:  Returns the median (unaligned) length of 
 *           the sequences in an alignment.
 */
static int
msa_median_length(ESL_MSA *msa)
{
  int  status;
  int *len;
  int  i;
  int  median;
  ESL_SQ *sq;
  sq = esl_sq_CreateDigital(msa->abc);

  ESL_ALLOC(len, sizeof(int) * msa->nseq);
  for (i = 0; i < msa->nseq; i++) {
    esl_sq_GetFromMSA(msa, i, sq);
    len[i] = sq->n;
    esl_sq_Reuse(sq);
    /*printf("i: %d len: %d\n", i, len[i]);*/
  }

  qsort(len, msa->nseq, sizeof(int), compare_ints);

  median = len[msa->nseq / 2];
  free(len);

  esl_sq_Destroy(sq);
  return median;

 ERROR:
  esl_fatal("msa_median_length() memory allocation error.");
  return 0.; /* NEVERREACHED */
}

/* Function: msa_remove_seqs_below_minlen()
 * 
 * Purpose:  Remove sequences in MSA whose dealigned length is less than a minimum length.
 */
static int
msa_remove_seqs_below_minlen(ESL_MSA *msa, float minlen, ESL_MSA **ret_new_msa)
{
  int  status;
  int *useme;
  int  i;

  ESL_MSA *new_msa;
  ESL_SQ *sq;
  sq = esl_sq_CreateDigital(msa->abc);

  ESL_ALLOC(useme, sizeof(int) * msa->nseq);
  for (i = 0; i < msa->nseq; i++) {
    esl_sq_GetFromMSA(msa, i, sq);
    useme[i] = ((float) sq->n >= minlen) ? TRUE : FALSE;
    /*printf("useme[i:%d]: %d\n", i, useme[i]);*/
    esl_sq_Reuse(sq);
  }

  if((status = esl_msa_SequenceSubset(msa, useme, &new_msa)) != eslOK) esl_fatal("esl_msa_SequenceSubset() had a problem.");
  free(useme);
  esl_sq_Destroy(sq);
  *ret_new_msa = new_msa;
  return eslOK;

 ERROR:
  esl_fatal("msa_remove_seqs_below_minlen() memory allocation error.");
  return eslOK; /* NEVERREACHED */
}

/* Function: msa_remove_truncated_seqs()
 * 
 * Purpose:  Remove sequences in MSA that have all gaps in the first <ntrunc> 5' leading 
 *           non-gap RF columns OR the last <ntrunc> 3' leading non-gap RF columns
 */
static int
msa_remove_truncated_seqs(ESL_MSA *msa, char *errbuf, int ntrunc, ESL_MSA **ret_new_msa)
{
  int  status;
  int *useme;
  int  i;
  int  leading_okay, trailing_okay;
  int  apos, cpos_ct;
  int  nused = 0;
  ESL_MSA *new_msa;

  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL)) ESL_XFAIL(eslEINVAL, errbuf, "in msa_remove_truncated_seqs(), msa must be digitized.");
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --detrunc.");

  ESL_ALLOC(useme, sizeof(int) * msa->nseq);

  for(i = 0; i < msa->nseq; i++) { 
    /* if ALL of the first 5' <ntrunc> non-gap RF columns are gaps in this seq, we'll remove it */
    leading_okay  = FALSE;
    cpos_ct = 0; 
    apos = 1;
    while(!leading_okay && (cpos_ct < ntrunc) && (apos <= msa->alen)) { 
      if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) { 
	cpos_ct++;
	if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) leading_okay = TRUE;
      }
      apos++;
    }

    trailing_okay = FALSE;
    cpos_ct = 0;
    apos = msa->alen;
    while(!trailing_okay && (cpos_ct < ntrunc) && (apos >= 1)) { 
      if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) { 
	cpos_ct++;
	if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) trailing_okay = TRUE;
      }
      apos--;
    }
    useme[i] = (leading_okay && trailing_okay) ? TRUE : FALSE;
    if(useme[i]) nused++;
  }
  if(nused == 0) ESL_FAIL(eslEINVAL, errbuf, "--detrunc removed ALL sequences!");
  if((status = esl_msa_SequenceSubset(msa, useme, &new_msa)) != eslOK) esl_fatal("esl_msa_SequenceSubset() had a problem.");
  free(useme);
  *ret_new_msa = new_msa;
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "msa_remove_truncated_seqs(): memory allocation error.");
  return eslOK; /* NEVERREACHED */
}

/* dump_infocontent
 *                   
 * Given an MSA with RF annotation, print a postscript heatmap of the 
 * information content of each non-gap RF column (consensus column).
 */
static int dump_infocontent(FILE *fp, ESL_MSA *msa, char *errbuf)
{
  int status;
  int apos, cpos;
  int i;
  int clen;
  double *obs, *ent, *bg;

  /* contract check */
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --icinfo.");
  if(! (msa->flags & eslMSA_DIGITAL))
    ESL_XFAIL(eslEINVAL, errbuf, "in dump_infocontent(), msa must be digitized.");

  clen = 0;
  for(apos = 1; apos <= msa->alen; apos++)
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) clen++;

  ESL_ALLOC(ent, sizeof(double) * clen);
  ESL_ALLOC(obs, sizeof(double) * msa->abc->K);
  ESL_ALLOC(bg, sizeof(double) * msa->abc->K);
  esl_vec_DSet(bg, msa->abc->K, 1./(msa->abc->K));

  cpos = 0;
  fprintf(fp, "# %4s  %5s\n", "cpos", "info");
  fprintf(fp, "# %4s  %5s\n", "----", "-----");
  for(apos = 1; apos <= msa->alen; apos++) {
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) { /* a consensus position */
      esl_vec_DSet(obs, msa->abc->K, 0.);
      for(i = 0; i < msa->nseq; i++)
	if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) { 
	  esl_abc_DCount(msa->abc, obs, msa->ax[i][apos], 1.);
	}
      esl_vec_DNorm(obs, msa->abc->K);
      ent[cpos] = esl_vec_DEntropy(bg, msa->abc->K) - esl_vec_DEntropy(obs, msa->abc->K);
      fprintf(fp, " %4d  %5.3f\n", cpos, ent[cpos]);
      cpos++;
    }
  }

  free(ent);
  free(obs);
  free(bg);
  return eslOK;

 ERROR:
  return status;
}

/* number_columns
 *                   
 * Add annotation to an MSA numbering the columns, either all
 * the columns (if <do_all>) or just non-gap #=GC RF columns.
 */
static int
number_columns(ESL_MSA *msa, int do_all, char *errbuf)
{
  int  status;
  int i;
  char *numstring;
  char *tag;
  int alen_ndigits;
  int tagwidth;
  int a,b,apos;
  int bmin;
  int pos2print;
  int tagidx;

  /* contract check */
  if(!do_all && msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment.");

  alen_ndigits = int_ndigits(msa->alen);
  tagwidth = do_all ? (3+alen_ndigits) : (5+alen_ndigits); /* "COL.X" or RFCOL.X" */

  ESL_ALLOC(tag, sizeof(char) * (tagwidth+1));
  ESL_ALLOC(numstring, sizeof(char) * (msa->alen+1));
  numstring[msa->alen] = '\0';
  tag[tagwidth] = '\0';
  if(do_all) { 
    bmin = 3;
    tag[0] = 'C';
    tag[1] = 'O';
    tag[2] = 'L';
  }
  else { 
    bmin = 5;
    tag[0] = 'R';
    tag[1] = 'F';
    tag[2] = 'C';
    tag[3] = 'O';
    tag[4] = 'L';
  }

  for(a = 0; a < alen_ndigits; a++) { 
    for(b = 0; b < alen_ndigits; b++) tag[b+bmin] = (a == b) ? 'X' : '.';
    pos2print = 1;
    for(apos = 1; apos <= msa->alen; apos++) { 
      if(!do_all && (esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)]))) numstring[(apos-1)] = '.';
      else numstring[(apos-1)] = get_char_digit_x_from_int(pos2print++, (alen_ndigits-a));
	/*printf("called get_char_digit_x_from_int(%d, %d)\n",apos, (alen_ndigits-a));*/
    }
    /* If the tag already exists, free it's associated markup string. This is an awful hack. */
    for (tagidx = 0; tagidx < msa->ngc; tagidx++) 
      if (strcmp(msa->gc_tag[tagidx], tag) == 0) break;
    if(tagidx != msa->ngc) { /* tag exists */
      free(msa->gc[tagidx]);
      msa->gc[tagidx] = NULL;
    }

    esl_msa_AppendGC(msa, tag, numstring);
  }

  ESL_ALLOC(numstring, sizeof(char) * (msa->alen + 1));
  for(i = 0; i < msa->alen; i++) { 
    numstring[i] = digit_to_char(i);
  }
  numstring[msa->alen] = '\0';
  free(numstring);
  return eslOK;

 ERROR:
  return eslEMEM;
}


/* digit_to_char
 *                   
 * Given a digit (0-9) return the character reprentation of it.
 * There must be a better way to do this; oh well.
 */
static char
digit_to_char(int digit) 
{
  if(digit == 0) return '0';
  if(digit == 1) return '1';
  if(digit == 2) return '2';
  if(digit == 3) return '3';
  if(digit == 4) return '4';
  if(digit == 5) return '5';
  if(digit == 6) return '6';
  if(digit == 7) return '7';
  if(digit == 8) return '8';
  if(digit == 9) return '9';
  else return '?';
}

/* Function: int_ndigits
 * Returns: The number of digits in <i>.
 */
static int
int_ndigits(int i)
{
  int n   = 0;
  while(i > 0) { i/=10; n++; }
  return n;
}

/* get_char_digit_x_from_int
 *                   
 * Given two integers <i> and <place> return the 
 * character version of the <place>'th digit in <i>.
 * Example <i> = 14378 <place> = 4 would return 7.
 */
static char
get_char_digit_x_from_int(int i, int place)
{
  int n,a,divisor;
  n = int_ndigits(i);

  if(n < place) return digit_to_char(0);

  divisor = 1;
  for(a = 0; a < (place-1); a++) divisor *= 10;
  /* subtract leading digits before the one we care about */
  i %= (divisor*10);
  return digit_to_char (i / divisor);
}

/* Function: read_seq_name_file
 * Date:     EPN, Thu Jun  5 13:21:36 2008
 * 
 * Read a file listing sequence names to remove or keep.
 * Store sequences in *ret_seqlist and return it.
 * Each white-space delimited token is considered a 
 * different sequence name. No checking is done in this 
 * function, but rather in subsequent functions. Each sequence name is 
 * 
 * Returns eslOK on success.
 */
int
read_seq_name_file(char *filename, char *errbuf, char ***ret_seqlist, int *ret_seqlist_n)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  int             toklen;
  int nalloc     = 10;
  int chunksize  = 10;
  char **seqlist = NULL;
  int n = 0;
  int i;
  void *tmp;

  ESL_ALLOC(seqlist, sizeof(char *) * nalloc);
  if (esl_fileparser_Open(filename, &efp) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "failed to open %s in read_seq_name_file\n", filename);
  
  while((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslEOF) {
    if(n == nalloc) { nalloc += chunksize; ESL_RALLOC(seqlist, tmp, sizeof(char *) * nalloc); }
    if((status = esl_strdup(tok, -1, &(seqlist[n++]))) != eslOK) ESL_FAIL(status, errbuf, "error in esl_strdup.");
  }
  esl_fileparser_Close(efp);
  *ret_seqlist = seqlist;
  *ret_seqlist_n = n;
  return eslOK;

 ERROR:
  if(seqlist != NULL) {
    for(i = 0; i < n; i++) free(seqlist[i]); 
    free(seqlist);
  }
  return status;
}


/* Function: msa_keep_or_remove_seqs()
 * 
 * Purpose:  Given a list of <seqlist_n> sequences in <seqlist>, either remove those
 *           sequences from msa, or remove all other sequences besides those from msa.
 *           Create and return the new msa with only the specified seqs in <ret_new_msa>.
 * 
 * Note:     Terribly inefficient, does a linear search for each seq, no sorting or anything.
 *
 * Returns: eslOK on success, eslEINVAL if a sequence name in seqlist does not exist in the msa.
 * 
 */
static int
msa_keep_or_remove_seqs(ESL_MSA *msa, char *errbuf, char **seqlist, int seqlist_n, int do_keep, ESL_MSA **ret_new_msa)
{
  int  status;
  int *useme;
  int  i, n;

  ESL_MSA *new_msa;
  ESL_SQ *sq;
  sq = esl_sq_CreateDigital(msa->abc);

  ESL_ALLOC(useme, sizeof(int) * msa->nseq);
  if(do_keep) esl_vec_ISet(useme, msa->nseq, FALSE);
  else        esl_vec_ISet(useme, msa->nseq, TRUE); 
  for(n = 0; n < seqlist_n; n++) { 
    for (i = 0; i < msa->nseq; i++) {
      if(strcmp(seqlist[n], msa->sqname[i]) == 0) { 
	useme[i] = do_keep ? TRUE : FALSE;
	break;
      }
      if(i == (msa->nseq-1)) ESL_FAIL(eslEINVAL, errbuf, "ERROR sequence %s does not exist in the MSA!", seqlist[n]);
    }
  }      

  if((status = esl_msa_SequenceSubset(msa, useme, &new_msa)) != eslOK) esl_fatal("esl_msa_SequenceSubset() had a problem.");
  free(useme);
  *ret_new_msa = new_msa;
  return eslOK;

 ERROR:
  esl_fatal("msa_keep_or_remove_seqs() memory allocation error.");
  return eslOK; /* NEVERREACHED */
}


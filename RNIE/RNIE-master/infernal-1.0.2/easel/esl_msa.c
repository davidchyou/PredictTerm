/*::cexcerpt::header_example::begin::*/
/* Multiple sequence alignment file i/o.
 *    
 * Contents:   
 *    1. The <ESL_MSA> object
 *    2. The <ESL_MSAFILE> object
 *    3. Digital mode MSA's         (augmentation: alphabet)
 *    4. Random MSA database access (augmentation: ssi)
 *    5. General i/o API, for all alignment formats
 *    6. Miscellaneous functions for manipulating MSAs
 *    7. Debugging/development routines
 *    8. Benchmark driver
 *    9. Unit tests
 *   10. Test driver
 *   11. Examples
 *   12. Copyright and license information
 *   
 * Augmentations:
 *   alphabet:  adds support for digital MSAs
 *   keyhash:   speeds up Stockholm file input
 *   ssi:       enables indexed random access in a file of many MSAs
 * 
 * SRE, Thu Jan 20 08:50:43 2005 [St. Louis]
 * SVN $Id: esl_msa.c 404 2009-10-12 13:54:50Z nawrockie $
 */
/*::cexcerpt::header_example::end::*/

/*::cexcerpt::include_example::begin::*/
#include "esl_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "easel.h"
#ifdef eslAUGMENT_KEYHASH
#include "esl_keyhash.h"
#endif
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#ifdef eslAUGMENT_SSI
#include "esl_ssi.h"
#endif
#include "esl_msa.h"
#include "esl_wuss.h"
/*::cexcerpt::include_example::end::*/



/******************************************************************************
 *# 1. The <ESL_MSA> object                                           
 *****************************************************************************/

/* create_mostly()
 * SRE, Sun Aug 27 16:40:00 2006 [Leesburg]
 *
 * This is the routine called by esl_msa_Create() and esl_msa_CreateDigital()
 * that does all allocation except the aseq/ax alignment data.
 * 
 * <nseq> may be the exact known # of seqs in an alignment; or <nseq>
 * may be an allocation block size (to be expanded by doubling, in
 * esl_msa_Expand(), as in:
 *     <if (msa->nseq == msa->sqalloc) esl_msa_Expand(msa);>
 * 
 * <alen> may be the exact length of an alignment, in columns; or it
 * may be -1, which states that your parser will take responsibility
 * for expanding as needed as new input is read into a growing new
 * alignment.
 *
 * A created <msa> can only be <_Expand()>'ed if <alen> is -1.
 *
 * Args:     <nseq> - number of sequences, or nseq allocation blocksize
 *           <alen> - length of alignment in columns, or -1     
 *
 * Returns:   pointer to new MSA object, w/ all values initialized.
 *            Note that msa->nseq is initialized to 0 here, even though space
 *            is allocated.
 *           
 * Throws:    <NULL> on allocation failure.          
 */
static ESL_MSA *
create_mostly(int nseq, int64_t alen)
{
  int      status;
  ESL_MSA *msa     = NULL;
  int      i;

  ESL_ALLOC(msa, sizeof(ESL_MSA));
  msa->aseq    = NULL;
  msa->sqname  = NULL;
  msa->wgt     = NULL;
  msa->alen    = alen;		/* if -1, then we're growable. */
  msa->nseq    = 0;		/* our caller (text or digital allocation) sets this.  */
  msa->flags   = 0;

#ifdef eslAUGMENT_ALPHABET
  msa->abc     = NULL;
  msa->ax      = NULL;
#endif /*eslAUGMENT_ALPHABET*/

  msa->name    = NULL;
  msa->desc    = NULL;
  msa->acc     = NULL;
  msa->au      = NULL;
  msa->ss_cons = NULL;
  msa->sa_cons = NULL;
  msa->rf      = NULL;
  msa->sqacc   = NULL;
  msa->sqdesc  = NULL;
  msa->ss      = NULL;
  msa->sa      = NULL;
  for (i = 0; i < eslMSA_NCUTS; i++) {
    msa->cutoff[i] = 0.;
    msa->cutset[i] = FALSE;
  }
  msa->sqalloc = nseq;
  msa->sqlen   = NULL;
  msa->sslen   = NULL;
  msa->salen   = NULL;
  msa->lastidx = 0;

  /* Unparsed markup, including comments and Stockholm tags.
   * GS, GC, and GR Stockholm tags require keyhash augmentation
   */
  msa->comment        = NULL;
  msa->ncomment       = 0;
  msa->alloc_ncomment = 0;

  msa->gf_tag         = NULL;
  msa->gf             = NULL;
  msa->ngf            = 0;
  msa->alloc_ngf      = 0;

  msa->gs_tag         = NULL;
  msa->gs             = NULL;
  msa->ngs            = 0;

  msa->gc_tag         = NULL;
  msa->gc             = NULL;
  msa->ngc            = 0;

  msa->gr_tag         = NULL;
  msa->gr             = NULL;
  msa->ngr            = 0;

#ifdef eslAUGMENT_KEYHASH
  msa->index     = esl_keyhash_Create();
  msa->gs_idx    = NULL;
  msa->gc_idx    = NULL;
  msa->gr_idx    = NULL;
#endif /*eslAUGMENT_KEYHASH*/

#ifdef eslAUGMENT_SSI
  msa->offset    = 0;
#endif

  /* Allocation, round 2.
   */
  ESL_ALLOC(msa->sqname, sizeof(char *) * nseq);
  ESL_ALLOC(msa->wgt,    sizeof(double) * nseq);
  ESL_ALLOC(msa->sqlen,  sizeof(int64_t)* nseq);

  /* Initialize at the second level.
   */
  for (i = 0; i < nseq; i++)
    {
      msa->sqname[i] = NULL;
      msa->sqlen[i]  = 0;
      msa->wgt[i]    = -1.0;	/* "unset so far" */
    }

  return msa;

 ERROR:
  esl_msa_Destroy(msa);
  return NULL;
}

/* get_seqidx()
 * 
 * Find the index of a given sequence <name> in an <msa>.
 * If caller has a good guess (for instance, the sequences are
 * coming in a previously seen order in a block of seqs or annotation),
 * the caller can pass this information in <guess>, or -1 if
 * it has no guess.
 * 
 * This function behaves differently depending on whether
 * keyhash augmentation is available or not. Without keyhashing,
 * the name is identified by bruteforce search of the names
 * in the <msa>. With keyhashing, we hash search, which should
 * improve performance for large alignments.
 * 
 * If the name does not already exist in the MSA, then it
 * is assumed to be a new sequence name that we need to store.
 * seqidx is set to msa->nseq, the MSA is Expand()'ed if necessary
 * to make room, the name is stored in msa->sqname[msa->nseq],
 * (and in the hash table, if we're keyhash augmented)
 * and msa->nseq is incremented.
 *
 * Returns:  <eslOK> on success, and the seqidx is 
 *           passed back via <ret_idx>. If <name> is new
 *           in the <msa>, the <name> is stored and the <msa> 
 *           may be internally reallocated if needed.
 *           
 * Throws:   <eslEMEM> if we try to add a name and allocation fails.
 *           <eslEINVAL> if we try to add a name to a non-growable MSA.
 */
static int
get_seqidx(ESL_MSA *msa, char *name, int guess, int *ret_idx)
{
  int seqidx;
  int status;

  *ret_idx = -1;

  /* can we guess? */
  if (guess >= 0 && 
      guess < msa->nseq && 
      strcmp(name, msa->sqname[guess]) == 0) 
    { *ret_idx = guess; return eslOK; }

  /* Else look it up - either brute force
   * or, if we're keyhash-augmented, by hashing.
   */
#ifdef eslAUGMENT_KEYHASH                  
  status = esl_key_Store(msa->index, name, &seqidx);
  if (status == eslEDUP) { *ret_idx = seqidx; return eslOK; }
  if (status != eslOK) return status; /* an error. */
#else
  for (seqidx = 0; seqidx < msa->nseq; seqidx++)
    if (strcmp(msa->sqname[seqidx], name) == 0) break;
  if (seqidx < msa->nseq) 
    { *ret_idx = seqidx; return eslOK; }
#endif

  /* If we reach here, then this is a new name that we're
   * adding.
   */
  if (seqidx >= msa->sqalloc &&  
     (status = esl_msa_Expand(msa)) != eslOK)
    return status; 
    
  status = esl_strdup(name, -1, &(msa->sqname[seqidx]));
  msa->nseq++;
  if (ret_idx != NULL) *ret_idx = seqidx;
  return status;
}

/* set_seq_accession()
 *
 * Sets the sequence accession field for sequence
 * number <seqidx> in an alignment <msa>, by
 * duplicating the string <acc>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
static int
set_seq_accession(ESL_MSA *msa, int seqidx, const char *acc)
{
  int status;
  int i;

  /* If this is the first acccession, we have to
   * initialize the whole optional array.
   */
  if (msa->sqacc == NULL) 
    {
      ESL_ALLOC(msa->sqacc, sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++)
	msa->sqacc[i] = NULL;
    }
  /* If we already had an accession, that's weird, but free it. 
   */
  if (msa->sqacc[seqidx] != NULL) free(msa->sqacc[seqidx]);
  return (esl_strdup(acc, -1, &(msa->sqacc[seqidx])));

 ERROR:
  return status;
}

/* set_seq_description()
 *
 * Set the sequence description field for sequence number
 * <seqidx> in an alignment <msa> by copying the string <desc>.
 *
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.
 */
static int
set_seq_description(ESL_MSA *msa, int seqidx, const char *desc)
{
  int status;
  int i;

  if (msa->sqdesc == NULL) 
    {
      ESL_ALLOC(msa->sqdesc, sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++)
	msa->sqdesc[i] = NULL;
  }
  if (msa->sqdesc[seqidx] != NULL) free(msa->sqdesc[seqidx]);
  return (esl_strdup(desc, -1, &(msa->sqdesc[seqidx])));

 ERROR:
  return status;
}


/* set_seq_ss() 
 *
 * Set the secondary structure annotation for sequence number
 * <seqidx> in an alignment <msa> by copying the string <ss>.
 *
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.
 */
static int
set_seq_ss(ESL_MSA *msa, int seqidx, const char *ss)
{
  int status;
  int i;

  if (msa->ss == NULL) 
    {
      ESL_ALLOC(msa->ss, sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++)
	msa->ss[i] = NULL;
  }
  if (msa->ss[seqidx] != NULL) free(msa->ss[seqidx]);
  return (esl_strdup(ss, -1, &(msa->ss[seqidx])));

 ERROR:
  return status;
}

/* set_seq_sa() 
 *
 * Set the sruface accessibility annotation for sequence number
 * <seqidx> in an alignment <msa> by copying the string <sa>.
 *
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.
 */
static int
set_seq_sa(ESL_MSA *msa, int seqidx, const char *sa)
{
  int status;
  int i;

  if (msa->sa == NULL) 
    {
      ESL_ALLOC(msa->sa, sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++)
	msa->sa[i] = NULL;
  }
  if (msa->sa[seqidx] != NULL) free(msa->sa[seqidx]);
  return (esl_strdup(sa, -1, &(msa->sa[seqidx])));

 ERROR:
  return status;
}


/* verify_parse()
 *
 * Last function called after a multiple alignment parser thinks it's
 * done. Checks that parse was successful; makes sure required
 * information is present; makes sure required information is
 * consistent. Some fields that are only use during parsing may be
 * freed (sqlen, for example), and some fields are finalized now
 * (<msa->alen> is set, for example). 
 * 
 * <errbuf> is a place to sprintf an informative message about the
 * reason for a parse error. The caller provides an <errbuf>
 * of at least 512 bytes.
 *
 * Returns:  <eslOK>, and errbuf is set to an empty string.
 *           
 * Throws:   <eslEFORMAT> if a problem is detected, and an
 *           informative message about the failure is in errbuf.
 */
static int
verify_parse(ESL_MSA *msa, char *errbuf)
{
  int idx;

  if (msa->nseq == 0) 
    {
      sprintf(errbuf, 
	      "MSA parse error: no sequences were found for alignment %.128s",
	      msa->name != NULL ? msa->name : "");
      return eslEFORMAT;
    }

  /* set alen, until proven otherwise; we'll check that the other seqs
   * have the same length later.
   */
  msa->alen = msa->sqlen[0];

  /* We can rely on msa->sqname[] being valid for any index,
   * because of the way the line parsers always store any name
   * they add to the index.
   */
  for (idx = 0; idx < msa->nseq; idx++)
    {
#ifdef eslAUGMENT_ALPHABET
      if ((msa->flags & eslMSA_DIGITAL) &&  (msa->ax  == NULL || msa->ax[idx] == NULL))
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: no sequence for %s",
		 msa->name != NULL ? msa->name : "", msa->sqname[idx]); 
#endif
      if (! (msa->flags & eslMSA_DIGITAL) && (msa->aseq == NULL || msa->aseq[idx] == NULL))
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: no sequence for %s",
		 msa->name != NULL ? msa->name : "", msa->sqname[idx]); 

      /* either all weights must be set, or none of them */
      if ((msa->flags & eslMSA_HASWGTS) && msa->wgt[idx] == -1.0)
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: expected a weight for seq %s", 
		  msa->name != NULL ? msa->name : "", msa->sqname[idx]);

      /* all aseq must be same length. */
      if (msa->sqlen[idx] != msa->alen)
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: sequence %s: length %" PRId64 ", expected %" PRId64,
		 msa->name != NULL ? msa->name : "", msa->sqname[idx], msa->sqlen[idx], msa->alen);

      /* if individual SS is present, it must have length right too */
      if (msa->ss != NULL &&  msa->ss[idx] != NULL &&  msa->sslen[idx] != msa->alen) 
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GR SS for %s: length %" PRId64 ", expected %" PRId64,
		 msa->name != NULL ? msa->name : "", msa->sqname[idx], msa->sslen[idx], msa->alen);

				/* if SA is present, must have length right */
      if (msa->sa != NULL && msa->sa[idx] != NULL && msa->salen[idx] != msa->alen) 
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GR SA for %s: length %" PRId64 ", expected %" PRId64,
		 msa->name != NULL ? msa->name : "", msa->sqname[idx], msa->salen[idx], msa->alen);
    }

  /* if cons SS is present, must have length right */
  if (msa->ss_cons != NULL && strlen(msa->ss_cons) != msa->alen) 
    ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GC SS_cons markup: len %zd, expected %" PRId64,
	     msa->name != NULL ? msa->name : "",  strlen(msa->ss_cons), msa->alen);

  /* if cons SA is present, must have length right */
  if (msa->sa_cons != NULL && strlen(msa->sa_cons) != msa->alen) 
    ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GC SA_cons markup: len %zd, expected %" PRId64,
	     msa->name != NULL ? msa->name : "",  strlen(msa->sa_cons), msa->alen);

  /* if RF is present, must have length right */
  if (msa->rf != NULL && strlen(msa->rf) != msa->alen) 
    ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GC RF markup: len %zd, expected %" PRId64,
	     msa->name != NULL ? msa->name : "", strlen(msa->rf), msa->alen);

  /* If no weights were set, set 'em all to 1.0 */
  if (!(msa->flags & eslMSA_HASWGTS))
    for (idx = 0; idx < msa->nseq; idx++)
      msa->wgt[idx] = 1.0;

  /* Clean up a little from the parser */
  if (msa->sqlen != NULL) { free(msa->sqlen); msa->sqlen = NULL; }
  if (msa->sslen != NULL) { free(msa->sslen); msa->sslen = NULL; }
  if (msa->salen != NULL) { free(msa->salen); msa->salen = NULL; }
  return eslOK;
}


/* Function:  esl_msa_Create()
 * Synopsis:  Creates an <ESL_MSA> object.
 * Incept:    SRE, Sun Jan 23 08:25:26 2005 [St. Louis]
 *
 * Purpose:   Creates and initializes an <ESL_MSA> object, and returns a
 *            pointer to it. 
 *  
 *            If caller already knows the dimensions of the alignment,
 *            both <nseq> and <alen>, then <msa = esl_msa_Create(nseq,
 *            alen)> allocates the whole thing at once. The MSA's
 *            <nseq> and <alen> fields are set accordingly, and the
 *            caller doesn't have to worry about setting them; it can
 *            just fill in <aseq>.
 *            
 *            If caller doesn't know the dimensions of the alignment
 *            (for example, when parsing an alignment file), then
 *            <nseq> is taken to be an initial allocation size, and
 *            <alen> must be -1. <alen=-1> is used as a flag for a
 *            "growable" MSA. For example, the call <msa =
 *            esl_msa_Create(16, -1)>.  allocates internally for an
 *            initial block of 16 sequences, but without allocating
 *            any space for individual sequences.  This allocation can
 *            be expanded (by doubling) by calling <esl_msa_Expand()>.
 *            A created <msa> can only be <_Expand()>'ed if <alen> is
 *            -1.
 *            
 *            In a growable alignment, caller becomes responsible for
 *            memory allocation of each individual <aseq[i]>. Caller
 *            is also responsible for setting <nseq> and <alen> when
 *            it is done parsing and creating the new MSA. In
 *            particular, the <esl_msa_Destroy()> function relies on
 *            <nseq> to know how many individual sequences are
 *            allocated.
 *
 * Args:      <nseq> - number of sequences, or nseq allocation blocksize
 *            <alen> - length of alignment in columns, or -1      
 *
 * Returns:   pointer to new MSA object, w/ all values initialized.
 *           
 * Throws:    <NULL> on allocation failure.          
 *
 * Xref:      squid's MSAAlloc()
 */
ESL_MSA *
esl_msa_Create(int nseq, int64_t alen)
{
  int      status;
  ESL_MSA *msa;
  int      i;

  msa = create_mostly(nseq, alen); /* aseq is null upon successful return */
  if (msa == NULL) return NULL; /* already threw error in mostly_create, so percolate */

  ESL_ALLOC(msa->aseq,   sizeof(char *) * msa->sqalloc);
  for (i = 0; i < msa->sqalloc; i++)
    msa->aseq[i] = NULL;

  if (alen != -1) {
    for (i = 0; i < nseq; i++)
      {
	ESL_ALLOC(msa->aseq[i], sizeof(char) * (alen+1));
	msa->aseq[i][alen] = '\0'; /* caller might forget to null terminate; help the poor */
      }
    msa->nseq = nseq;
  }
  return msa;

 ERROR:
  esl_msa_Destroy(msa);
  return NULL;
}


/* Function:  esl_msa_Expand()
 * Synopsis:  Reallocate for more sequences.
 * Incept:    SRE, Sun Jan 23 08:26:30 2005 [St. Louis]
 *
 * Purpose:   Double the current sequence allocation in <msa>.
 *            Typically used when we're reading an alignment sequentially 
 *            from a file, so we don't know nseq 'til we're done.
 *            
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEMEM> on reallocation failure; <msa> is undamaged,
 *            and the caller may attempt to recover from the error.
 *            
 *            Throws <eslEINVAL> if <msa> is not growable: its <alen>
 *            field must be -1 to be growable.
 *
 * Xref:      squid's MSAExpand(), 1999.
 */
int
esl_msa_Expand(ESL_MSA *msa)
{
  int   status;
  int   old, new;		/* old & new allocation sizes (max # seqs) */
  void *p;			/* tmp ptr to realloc'ed memory */
  int   i,j;

  if (msa->alen != -1) 
    ESL_EXCEPTION(eslEINVAL, "that MSA is not growable");

  old = msa->sqalloc;
  new = 2*old;

  /* Normally either aseq (ascii) or ax (digitized) would be active, not both.
   * We could make sure that that's true, but that's checked elsewhere.           
   */
  if (msa->aseq != NULL) ESL_RALLOC(msa->aseq, p, sizeof(char *)    * new);
#ifdef eslAUGMENT_ALPHABET
  if (msa->ax   != NULL) ESL_RALLOC(msa->ax,   p, sizeof(ESL_DSQ *) * new);
#endif /*eslAUGMENT_ALPHABET*/

  ESL_RALLOC(msa->sqname, p, sizeof(char *) * new);
  ESL_RALLOC(msa->wgt,    p, sizeof(double) * new);
  ESL_RALLOC(msa->sqlen,  p, sizeof(int64_t)* new);

  if (msa->ss != NULL) 
    {
      ESL_RALLOC(msa->ss,    p, sizeof(char *)  * new);
      ESL_RALLOC(msa->sslen, p, sizeof(int64_t) * new);
    }
  
  if (msa->sa != NULL) 
    {
      ESL_RALLOC(msa->sa,    p, sizeof(char *)  * new);
      ESL_RALLOC(msa->salen, p, sizeof(int64_t) * new);
    }

  if (msa->sqacc != NULL)
    ESL_RALLOC(msa->sqacc,  p, sizeof(char *) * new);

  if (msa->sqdesc != NULL)
    ESL_RALLOC(msa->sqdesc, p, sizeof(char *) * new);

  for (i = old; i < new; i++)
    {
      if (msa->aseq != NULL) msa->aseq[i] = NULL;
#ifdef eslAUGMENT_ALPHABET
      if (msa->ax   != NULL) msa->ax[i]   = NULL;
#endif /*eslAUGMENT_ALPHABET*/
      msa->sqname[i] = NULL;
      msa->wgt[i]    = -1.0;	/* -1.0 means "unset so far" */
      msa->sqlen[i]  = 0;

      if (msa->ss != NULL) 
	{
	  msa->ss[i] = NULL;
	  msa->sslen[i] = 0;
	}
      if (msa->sa != NULL) 
	{ 
	  msa->sa[i] = NULL;
	  msa->salen[i] = 0;
	}
      if (msa->sqacc  != NULL) msa->sqacc[i]  = NULL;
      if (msa->sqdesc != NULL) msa->sqdesc[i] = NULL;
    }

  /* Reallocate and re-init for unparsed #=GS tags, if we have some.
   * gs is [0..ngs-1][0..nseq-1][], so we're reallocing the middle
   * set of pointers.
   */
  if (msa->gs != NULL)
    for (i = 0; i < msa->ngs; i++)
      {
	if (msa->gs[i] != NULL)
	  {
	    ESL_RALLOC(msa->gs[i], p, sizeof(char *) * new);
	    for (j = old; j < new; j++)
	      msa->gs[i][j] = NULL;
	  }
      }
  /* Reallocate and re-init for unparsed #=GR tags, if we have some.
   * gr is [0..ngs-1][0..nseq-1][], so we're reallocing the middle
   * set of pointers.
   */
  if (msa->gr != NULL)
    for (i = 0; i < msa->ngr; i++)
      {
	if (msa->gr[i] != NULL)
	  {
	    ESL_RALLOC(msa->gr[i], p, sizeof(char *) * new);
	    for (j = old; j < new; j++)
	      msa->gr[i][j] = NULL;
	  }
      }

  msa->sqalloc = new;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_msa_Copy()
 * Synopsis:  Copies an MSA.
 * Incept:    SRE, Tue Jan 22 15:30:32 2008 [Janelia]
 *
 * Purpose:   Makes a copy of <msa> in <new>. Caller has
 *            already allocated <new> to hold an MSA of
 *            at least <msa->nseq> sequences and <msa->alen>
 *            columns.
 *            
 * Note:      Because MSA's are not reusable, this function does a
 *            lot of internal allocation for optional fields, without
 *            checking <new> to see if space was already allocated. To
 *            reuse an MSA <new> and copy new data into it, we'll
 *            eventually need a <esl_msa_Reuse()> function, and/or
 *            recode this to reuse or free any already-allocated
 *            optional memory it encounters in <new>. Until then, 
 *            it's unlikely that <esl_msa_Copy()> is useful on its own;
 *            the caller would be expected to call <esl_msa_Clone()> 
 *            instead.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. In this case, <new>
 *            was only partially constructed, and should be treated
 *            as corrupt.
 */
int
esl_msa_Copy(const ESL_MSA *msa, ESL_MSA *new)
{
  int i, x, j;
  int status;

  /* aseq[0..nseq-1][0..alen-1] strings,
   * or ax[0..nseq-1][(0) 1..alen (alen+1)] digital seqs 
   * <new> must have one of them allocated already.
   */
  if (! (msa->flags & eslMSA_DIGITAL))
    for (i = 0; i < msa->nseq; i++)
      strcpy(new->aseq[i], msa->aseq[i]);
#ifdef eslAUGMENT_ALPHABET
  else
    {
      for (i = 0; i < msa->nseq; i++)
	memcpy(new->ax[i], msa->ax[i], (msa->alen+2) * sizeof(ESL_DSQ));
      new->abc = msa->abc;
    }
#endif
  
  for (i = 0; i < msa->nseq; i++) {
    esl_strdup(msa->sqname[i], -1, &(new->sqname[i]));
    new->wgt[i] = msa->wgt[i];
  }
  /* alen, nseq were already set by Create() */
  new->flags = msa->flags;

  esl_strdup(msa->name, -1, &(new->name));
  esl_strdup(msa->desc, -1, &(new->desc));
  esl_strdup(msa->acc,  -1, &(new->acc));
  esl_strdup(msa->au,    -1, &(new->au));
  esl_strdup(msa->ss_cons,   -1, &(new->ss_cons));
  esl_strdup(msa->sa_cons,   -1, &(new->sa_cons));
  esl_strdup(msa->rf,    -1, &(new->rf));

  if (msa->sqacc != NULL) {
    ESL_ALLOC(new->sqacc, sizeof(char **) * msa->nseq);
    for (i = 0; i < msa->nseq; i++)
      esl_strdup(msa->sqacc[i], -1, &(new->sqacc[i]));
  }
  if (msa->sqdesc != NULL) {
    ESL_ALLOC(new->sqdesc, sizeof(char **) * msa->nseq);
    for (i = 0; i < msa->nseq; i++)
      esl_strdup(msa->sqdesc[i], -1, &(new->sqdesc[i]));
  }
  if (msa->ss != NULL) {
    ESL_ALLOC(new->ss, sizeof(char **) * msa->nseq);
    for (i = 0; i < msa->nseq; i++)
      esl_strdup(msa->ss[i], -1, &(new->ss[i]));
  }
  if (msa->sa != NULL) {
    ESL_ALLOC(new->sa, sizeof(char **) * msa->nseq);
    for (i = 0; i < msa->nseq; i++)
      esl_strdup(msa->sa[i], -1, &(new->sa[i]));
  }
  
  for (x = 0; x < eslMSA_NCUTS; x++) {
    new->cutoff[x] = msa->cutoff[x];
    new->cutset[x] = msa->cutset[x];
  }

  if (msa->ncomment > 0) {
    ESL_ALLOC(new->comment, sizeof(char **) * msa->ncomment);
    new->ncomment       = msa->ncomment;
    new->alloc_ncomment = msa->ncomment;
    for (i = 0; i < msa->ncomment; i++)
      esl_strdup(msa->comment[i], -1, &(new->comment[i]));
  }

  if (msa->ngf > 0) {
    ESL_ALLOC(new->gf_tag, sizeof(char **) * msa->ngf);
    ESL_ALLOC(new->gf,     sizeof(char **) * msa->ngf);
    new->ngf       = msa->ngf;
    new->alloc_ngf = msa->ngf;
    for (i = 0; i < msa->ngf; i++) {
      esl_strdup(msa->gf_tag[i], -1, &(new->gf_tag[i]));
      esl_strdup(msa->gf[i],     -1, &(new->gf[i]));
    }
  }

  if (msa->ngs > 0) {
    ESL_ALLOC(new->gs_tag, sizeof(char **)  * msa->ngs);
    ESL_ALLOC(new->gs,     sizeof(char ***) * msa->ngs);
    new->ngs       = msa->ngs;
    for (i = 0; i < msa->ngs; i++) {
      ESL_ALLOC(new->gs[i], sizeof(char **) * msa->nseq);
      esl_strdup(msa->gs_tag[i], -1, &(new->gs_tag[i]));
      for (j = 0; j < msa->nseq; j++)
	esl_strdup(msa->gs[i][j],  -1, &(new->gs[i][j]));
    }
  }

  if (msa->ngc > 0) {
    ESL_ALLOC(new->gc_tag, sizeof(char **) * msa->ngc);
    ESL_ALLOC(new->gc,     sizeof(char **) * msa->ngc);
    new->ngc       = msa->ngc;
    for (i = 0; i < msa->ngc; i++) {
      esl_strdup(msa->gc_tag[i], -1, &(new->gc_tag[i]));
      esl_strdup(msa->gc[i],     -1, &(new->gc[i]));
    }
  }
  
  if (msa->ngr > 0) {
    ESL_ALLOC(new->gr_tag, sizeof(char **)  * msa->ngr);
    ESL_ALLOC(new->gr,     sizeof(char ***) * msa->ngr);
    new->ngr       = msa->ngr;
    for (i = 0; i < msa->ngr; i++) {
      ESL_ALLOC(new->gr[i], sizeof(char **) * msa->nseq);
      esl_strdup(msa->gr_tag[i], -1, &(new->gr_tag[i]));
      for (j = 0; j < msa->nseq; j++)
	esl_strdup(msa->gr[i][j],  -1, &(new->gr[i][j]));
    }
  }

#ifdef eslAUGMENT_KEYHASH
  if (msa->index  != NULL) new->index  = esl_keyhash_Clone(msa->index);
  if (msa->gs_idx != NULL) new->gs_idx = esl_keyhash_Clone(msa->gs_idx);
  if (msa->gc_idx != NULL) new->gc_idx = esl_keyhash_Clone(msa->gc_idx);
  if (msa->gr_idx != NULL) new->gr_idx = esl_keyhash_Clone(msa->gr_idx);
#endif

#ifdef eslAUGMENT_SSI
  new->offset = msa->offset;
#endif

  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_msa_Clone()
 * Synopsis:  Duplicates an MSA.
 * Incept:    SRE, Tue Jan 22 15:23:55 2008 [Janelia]
 *
 * Purpose:   Make a duplicate of <msa>, in newly 
 *            allocated space. 
 *
 * Returns:   a pointer to the newly allocated clone.
 *            Caller is responsible for free'ing it.
 *
 * Throws:    <NULL> on allocation error.
 */
ESL_MSA *
esl_msa_Clone(const ESL_MSA *msa)
{
  ESL_MSA *nw = NULL;
  int      status;

#ifdef eslAUGMENT_ALPHABET
  if (msa->flags & eslMSA_DIGITAL) {
      if ((nw = esl_msa_CreateDigital(msa->abc, msa->nseq, msa->alen)) == NULL)  return NULL;
  } else
#endif
  if ((nw     = esl_msa_Create(msa->nseq, msa->alen)) == NULL)  return NULL;  

  if ((status = esl_msa_Copy(msa, nw) )               != eslOK) goto ERROR;
  return nw;

 ERROR:
  esl_msa_Destroy(nw);
  return NULL;
}


/* Function:  esl_msa_Destroy()
 * Synopsis:  Frees an <ESL_MSA>.
 * Incept:    SRE, Sun Jan 23 08:26:02 2005 [St. Louis]
 *
 * Purpose:   Destroys <msa>.
 *
 * Xref:      squid's MSADestroy().
 */
void
esl_msa_Destroy(ESL_MSA *msa)
{
  if (msa == NULL) return;

  if (msa->aseq != NULL) 
    esl_Free2D((void **) msa->aseq, msa->nseq);
#ifdef eslAUGMENT_ALPHABET
  if (msa->ax != NULL) 
    esl_Free2D((void **) msa->ax, msa->nseq);
#endif /*eslAUGMENT_ALPHABET*/

  esl_Free2D((void **) msa->sqname, msa->nseq);
  esl_Free2D((void **) msa->sqacc,  msa->nseq);
  esl_Free2D((void **) msa->sqdesc, msa->nseq);
  esl_Free2D((void **) msa->ss,     msa->nseq);
  esl_Free2D((void **) msa->sa,     msa->nseq);

  if (msa->sqlen   != NULL) free(msa->sqlen);
  if (msa->wgt     != NULL) free(msa->wgt);

  if (msa->name    != NULL) free(msa->name);
  if (msa->desc    != NULL) free(msa->desc);
  if (msa->acc     != NULL) free(msa->acc);
  if (msa->au      != NULL) free(msa->au);
  if (msa->ss_cons != NULL) free(msa->ss_cons);
  if (msa->sa_cons != NULL) free(msa->sa_cons);
  if (msa->rf      != NULL) free(msa->rf);
  if (msa->sslen   != NULL) free(msa->sslen);
  if (msa->salen   != NULL) free(msa->salen);
  
  esl_Free2D((void **) msa->comment, msa->ncomment);
  esl_Free2D((void **) msa->gf_tag,  msa->ngf);
  esl_Free2D((void **) msa->gf,      msa->ngf);

  esl_Free2D((void **) msa->gs_tag,  msa->ngs);
  esl_Free3D((void ***)msa->gs,      msa->ngs, msa->nseq);
  esl_Free2D((void **) msa->gc_tag,  msa->ngc);
  esl_Free2D((void **) msa->gc,      msa->ngc);
  esl_Free2D((void **) msa->gr_tag,  msa->ngr);
  esl_Free3D((void ***)msa->gr,      msa->ngr, msa->nseq);

#ifdef eslAUGMENT_KEYHASH
  esl_keyhash_Destroy(msa->index);
  esl_keyhash_Destroy(msa->gs_idx);
  esl_keyhash_Destroy(msa->gc_idx);
  esl_keyhash_Destroy(msa->gr_idx);
#endif /* keyhash augmentation */  

  free(msa);
  return;
}


/* Function:  esl_msa_SetName()
 * Synopsis:  Format and set a name of an MSA.
 * Incept:    SRE, Sat Feb 23 18:42:47 2008 [Casa de Gatos]
 *
 * Purpose:   Sets the name of the msa <msa> to <name>, 
 *            (re-)allocating as needed.
 *
 *            <name> may be a <printf()>-style format with
 *            arguments; for example, <esl_msa_SetName(msa, "random%d", i)>.
 *
 *            A copy of <name> is made, so if caller 
 *            allocated <name>, it is still responsible 
 *            for free'ing it.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_msa_SetName(ESL_MSA *msa, const char *name, ...)
{
  va_list argp, argp2;
  int     n1 = 32;
  int     n2 = 0;
  void   *tmp;
  int     status;

  if (msa->name != NULL) free(msa->name);
  ESL_ALLOC(msa->name, sizeof(char) * n1);

  va_start(argp, name);
  va_copy(argp2, argp);
  if ((n2 = vsnprintf(msa->name, n1, name, argp)) >= n1)
    {
      ESL_RALLOC(msa->name, tmp, sizeof(char) * (n2+1));
      vsnprintf(msa->name, n2+1, name, argp2);
    }
  va_end(argp);
  va_end(argp2);
  return eslOK;

 ERROR:
  msa->name = NULL;
  return status;
}


/* Function:  esl_msa_SetDesc()
 * Synopsis:  Format and set description line of an MSA.
 * Incept:    SRE, Sat Feb 23 18:47:06 2008 [Casa de Gatos]
 *
 * Purpose:   Sets the description line of the msa <msa> to <desc>, 
 *            (re-)allocating as needed.
 *
 *            <desc> may be a <printf()>-style format with
 *            arguments; for example, <esl_msa_SetDesc(msa, "sample %d", i)>.
 *
 *            A copy of <desc> is made, so if caller 
 *            allocated <desc>, it is still responsible 
 *            for free'ing it.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_msa_SetDesc(ESL_MSA *msa, const char *desc, ...)
{
  va_list argp, argp2;
  int     n1 = 64;
  int     n2 = 0;
  void   *tmp;
  int     status;

  if (msa->desc != NULL) free(msa->desc);
  ESL_ALLOC(msa->desc, sizeof(char) * n1);

  va_start(argp, desc);
  va_copy(argp2, argp);
  if ((n2 = vsnprintf(msa->desc, n1, desc, argp)) >= n1)
    {
      ESL_RALLOC(msa->desc, tmp, sizeof(char) * (n2+1));
      vsnprintf(msa->desc, (n2+1), desc, argp2);
    }
  va_end(argp);
  va_end(argp2);
  return eslOK;

 ERROR:
  msa->desc = NULL;
  return status;
}

/* Function:  esl_msa_SetAccession()
 * Synopsis:  Format and set accession number of an MSA.
 * Incept:    SRE, Sat Feb 23 18:49:04 2008 [Casa de Gatos]
 *
 * Purpose:   Sets accession numbre of the msa <msa> to <acc>, 
 *            (re-)allocating as needed.
 *
 *            <acc> may be a <printf()>-style format with
 *            arguments; for example, <esl_msa_SetAccession(msa, "PF%06d", i)>.
 *
 *            A copy of <acc> is made, so if caller 
 *            allocated <acc>, it is still responsible 
 *            for free'ing it.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_msa_SetAccession(ESL_MSA *msa, const char *acc, ...)
{
  va_list argp, argp2;
  int     n1 = 16;
  int     n2 = 0;
  void   *tmp;
  int     status;

  if (msa->acc != NULL) free(msa->acc);
  ESL_ALLOC(msa->acc, sizeof(char) * n1);

  va_start(argp, acc);
  va_copy(argp2, argp);
  if ((n2 = vsnprintf(msa->acc, n1, acc, argp)) >= n1)
    {
      ESL_RALLOC(msa->acc, tmp, sizeof(char) * (n2+1));
      vsnprintf(msa->acc, (n2+1), acc, argp2);
    }
  va_end(argp);
  va_end(argp2);
  return eslOK;

 ERROR:
  msa->acc = NULL;
  return status;
}





/*---------------------- end of ESL_MSA functions ---------------------------*/



/******************************************************************************
 *# 2. The <ESL_MSAFILE> object                                       
 *****************************************************************************/

/* msafile_open():
 * this is the routine that actually opens an ESL_MSAFILE;
 * esl_msafile_Open() and esl_msafile_OpenDigital() are wrappers around it.
 */
static int
msafile_open(const char *filename, int format, const char *env, ESL_MSAFILE **ret_msafp)
{
  ESL_MSAFILE *afp     = NULL;
  char        *ssifile = NULL;
  char        *envfile = NULL;
  char        *cmd     = NULL;
  int          n       = strlen(filename);
  int          status;
  
  ESL_ALLOC(afp, sizeof(ESL_MSAFILE));
  afp->f          = NULL;
  afp->fname      = NULL;
  afp->linenumber = 0;
  afp->errbuf[0]  = '\0';
  afp->buf        = NULL;
  afp->buflen     = 0;
  afp->do_gzip    = FALSE;
  afp->do_stdin   = FALSE;
  afp->format     = eslMSAFILE_UNKNOWN;	
  afp->do_digital = FALSE;
#ifdef eslAUGMENT_ALPHABET
  afp->abc        = NULL;	        
#endif
#ifdef eslAUGMENT_SSI		
  afp->ssi        = NULL;	         
#endif  
  afp->msa_cache  = NULL;

  if (strcmp(filename, "-") == 0)
    {
      afp->f         = stdin;
      afp->do_stdin  = TRUE; 
      if ((status = esl_strdup("[STDIN]", -1, &(afp->fname))) != eslOK) goto ERROR;
    }
#ifdef HAVE_POPEN
  /* popen(), pclose() aren't portable to non-POSIX systems; 
   * disable this section in strict ANSI C mode.
   */
  /* tricky: if n= length of a string s, then
   * s+n-i repositions pointer s at the last i chars
   * of the string.
   */
  else if (n > 3 && strcmp(filename+n-3, ".gz") == 0)
    {
      /* Note that popen() will return "successfully"
       * if file doesn't exist, because gzip works fine
       * and prints an error! So we have to check for
       * existence of file ourself.
       */
      if (! esl_FileExists(filename))	      { status = eslENOTFOUND; goto ERROR; }
      ESL_ALLOC(cmd, sizeof(char) * (n+1+strlen("gzip -dc ")));
      sprintf(cmd, "gzip -dc %s", filename);
      if ((afp->f = popen(cmd, "r")) == NULL) { status = eslENOTFOUND; goto ERROR; }
      if ((status = esl_strdup(filename, n, &(afp->fname))) != eslOK)  goto ERROR;
      afp->do_gzip  = TRUE;
    }
#endif /*HAVE_POPEN*/
  else	/* Normal file open or env file open: set ssifile */
    {
      /* When we open a file, it may be either in the current
       * directory, or in the directory indicated by the env
       * argument - and we construct an SSI filename accordingly.
       * (Whether or not we're SSI augmented, in fact, for simplicity.)
       */
      if ((afp->f = fopen(filename, "r")) != NULL)
	{
	  if (esl_strdup(filename, n, &ssifile)                      != eslOK) goto ERROR;
	  if (esl_strcat(&ssifile, n, ".ssi", 4)                     != eslOK) goto ERROR;
	  if ((status = esl_strdup(filename, n, &(afp->fname)))       != eslOK) goto ERROR;
	}
      else if (esl_FileEnvOpen(filename, env, &(afp->f), &envfile) == eslOK)
	{
	  if (esl_strdup(envfile, n, &ssifile)                      != eslOK) goto ERROR;
	  if (esl_strcat(&ssifile, n, ".ssi", 4)                    != eslOK) goto ERROR;
	  if ((status = esl_strdup(envfile, n, &(afp->fname)))       != eslOK) goto ERROR;
	}
      else 
	{ status = eslENOTFOUND; goto ERROR;}
    }

#ifdef eslAUGMENT_SSI
  /* If augmented by SSI indexing:
   * Open the SSI index file. If it doesn't exist, or
   * it's corrupt, or some error happens, afp->ssi stays NULL.
   * We should warn, probably, or provide some way for caller to 
   * to know that we've opened the index successfully or not.
   */
  status = esl_ssi_Open(ssifile, &(afp->ssi)); 
#endif

  /* Invoke autodetection if we haven't already been told what
   * to expect.
   */
  if (format == eslMSAFILE_UNKNOWN)
    {
      if (afp->do_stdin == TRUE || afp->do_gzip)
	ESL_XEXCEPTION(eslEINVAL, "Can't autodetect alignment file fmt in stdin, gzip pipe");
      if (esl_msa_GuessFileFormat(afp) != eslOK)
	{ status = eslEFORMAT; goto ERROR; }
    }
  else 
    afp->format     = format;

  if (envfile != NULL) free(envfile);
  if (ssifile != NULL) free(ssifile);
  if (cmd     != NULL) free(cmd);
  *ret_msafp = afp;
  return eslOK;

 ERROR:
  if (envfile != NULL) free(envfile);
  if (ssifile != NULL) free(ssifile);
  if (cmd     != NULL) free(cmd);
  if (afp     != NULL) esl_msafile_Close(afp); 
  *ret_msafp = NULL;
  return status;
}


/* Function: esl_msafile_Open()
 * Synopsis: Open an MSA file for input.
 * Date:     SRE, Sun Jan 23 08:30:33 2005 [St. Louis]
 *
 * Purpose:  Open an alignment database file <filename> and prepare for
 *           reading one alignment, or sequentially in the case of 
 *           multiple MSA databases (e.g. Stockholm format); returns
 *           the opened file pointer in <ret_msafp>.
 *          
 *           There are one or two special cases for <filename>. If
 *           <filename> is "-", then the alignment is read from
 *           <stdin>. If <filename> ends in ".gz", then the file is
 *           assumed to be compressed by gzip, and it is opened as a
 *           pipe from <gunzip -dc>. (Auto-decompression of gzipp'ed files
 *           is only available on POSIX-compliant systems w/ popen(), when 
 *           <HAVE_POPEN> is defined at compile-time.)
 *          
 *           If <env> is non-NULL, then we look for <filename> in
 *           one or more directories in a colon-delimited list
 *           that is the value of the environment variable <env>.
 *           For example, if we had 
 *              <setenv HMMERDB /nfs/db/Pfam:/nfs/db/Rfam> 
 *           in the environment, a profile HMM application
 *           might pass "HMMERDB" as <env>.
 *          
 *          The file is asserted to be in format <fmt>, which is
 *          either a known format like <eslMSAFILE_STOCKHOLM>, or
 *          <eslMSAFILE_UNKNOWN>; if <fmt> is <eslMSAFILE_UNKNOWN>,
 *          then format autodetection is invoked.
 *
 * Returns:  <eslOK> on success, and <ret_msafp> is set to point at
 *           an open <ESL_MSAFILE>. Caller frees this file pointer with
 *           <esl_msafile_Close()>.
 *           
 *           Returns <eslENOTFOUND> if <filename> cannot be opened,
 *           or <eslEFORMAT> if autodetection is attempted and 
 *           format cannot be determined.
 *           
 * Throws:   <eslEMEM> on allocation failure.
 *           <eslEINVAL> if format autodetection is attempted on 
 *           stdin or a gunzip pipe.
 *
 * Xref:     squid's MSAFileOpen(), 1999.
 * 
 * Note      Implemented as a wrapper around msafile_open(), because
 *           esl_msafile_OpenDigital() shares almost all the same code.
 */
int
esl_msafile_Open(const char *filename, int format, const char *env, ESL_MSAFILE **ret_msafp)
{
  return msafile_open(filename, format, env, ret_msafp);
}

/* Function:  esl_msafile_Close()
 * Synopsis:  Closes an open MSA file.
 * Incept:    SRE, Sun Jan 23 08:18:39 2005 [St. Louis]
 *
 * Purpose:   Close an open <ESL_MSAFILE>.
 *
 * Xref:      squid's MSAFileClose().
 */
void
esl_msafile_Close(ESL_MSAFILE *afp)
{
  if (afp == NULL) return;

#ifdef HAVE_POPEN /* gzip functionality */
  if (afp->do_gzip && afp->f != NULL)    pclose(afp->f);
#endif
  if (!afp->do_gzip && ! afp->do_stdin && afp->f != NULL) fclose(afp->f);
  if (afp->fname != NULL) free(afp->fname);
  if (afp->buf  != NULL)  free(afp->buf);
#ifdef eslAUGMENT_SSI
  if (afp->ssi  != NULL)  esl_ssi_Close(afp->ssi); 
#endif /* eslAUGMENT_SSI*/
  if (afp->msa_cache != NULL) esl_msa_Destroy(afp->msa_cache);
  free(afp);
}

/*-------------------- end of ESL_MSAFILE functions -------------------------*/




/******************************************************************************
 *# 3. Digital mode MSA's (augmentation: alphabet)
 *****************************************************************************/
#ifdef eslAUGMENT_ALPHABET
/* Function:  esl_msa_GuessAlphabet()
 * Synopsis:  Guess alphabet of MSA.
 * Incept:    SRE, Fri May 18 09:55:08 2007 [Janelia]
 *
 * Purpose:   Guess whether the sequences in the <msa> are
 *            <eslDNA>, <eslRNA>, or <eslAMINO>, and return
 *            that guess in <*ret_type>.
 *            
 *            The determination is made based on the classifications
 *            of the individual sequences in the alignment. At least
 *            one sequence must contain ten residues or more to be
 *            classified. If one or more sequences is called
 *            <eslAMINO> and one or more is called <eslDNA>/<eslRNA>,
 *            the alignment's alphabet is considered to be
 *            indeterminate (<eslUNKNOWN>). If some sequences are
 *            <eslDNA> and some are <eslRNA>, the alignment is called
 *            <eslDNA>; this should cause no problems, because Easel
 *            reads U as a synonym for T in DNA sequence anyway.
 *            
 *            Tested on Pfam 21.0 and Rfam 7.0, this routine correctly
 *            classified all 8957 Pfam alignments as protein, and 503
 *            Rfam alignments as RNA (both seed and full alignments).
 *
 * Returns:   <eslOK> on success, and <*ret_type> is set
 *            to <eslDNA>, <eslRNA>, or <eslAMINO>. 
 *            
 *            Returns <eslEAMBIGUOUS> and sets <*ret_type> to
 *            <eslUNKNOWN> if the alphabet cannot be reliably guessed.
 *
 * Xref:      J1/62
 */
int
esl_msa_GuessAlphabet(const ESL_MSA *msa, int *ret_type)
{
  int64_t namino   = 0,
          ndna     = 0,
          nrna     = 0,
          nunknown = 0;
  int     type;
  int     i,x;
  int64_t j,n;
  int64_t ct[26];

  if (msa->flags & eslMSA_DIGITAL) { *ret_type = msa->abc->type; return eslOK; }

  *ret_type = eslUNKNOWN;

  /* On wide alignments, we're better off looking at individual sequence
   * classifications. We don't want to end up calling the whole alignment
   * indeterminate just because a few sequences have degenerate residue
   * codes.
   */
  for (i = 0; i < msa->nseq; i++) 
    {
      for (x = 0; x < 26; x++) ct[x] = 0;
      for (n = 0, j = 0; j < msa->alen; j++) {
	x = toupper(msa->aseq[i][j]) - 'A';
	if (x < 0 || x > 26) continue;
	ct[x]++;
	n++;
	if (n > 10000) break;	/* ought to know by now */
      }
      esl_abc_GuessAlphabet(ct, &type);

      switch (type) {
      case eslAMINO:   namino++; break;
      case eslDNA:     ndna++;   break;
      case eslRNA:     nrna++;   break;
      default:         nunknown++; 
      }
    }
  if      (namino    > 0 && (ndna+nrna)   == 0) *ret_type = eslAMINO;
  else if (ndna      > 0 && (nrna+namino) == 0) *ret_type = eslDNA;
  else if (nrna      > 0 && (ndna+namino) == 0) *ret_type = eslRNA;
  else if (ndna+nrna > 0 && namino        == 0) *ret_type = eslDNA;

  /* On narrow alignments, no single sequence may be long enough to 
   * be classified, but we can determine alphabet from composition
   * of the complete alignment. Of course, degenerate residue codes in
   * a DNA alignment will still screw us.
   */
  if (*ret_type == eslUNKNOWN)
    {

      n = 0;
      for (x = 0; x < 26; x++) ct[x] = 0;
      for (i = 0; i < msa->nseq; i++) {
	for (j = 0; j < msa->alen; j++) {
	  x = toupper(msa->aseq[i][j]) - 'A';
	  if (x < 0 || x > 26) continue;
	  ct[x]++;
	  n++;
	  if (n > 10000) break;	/* ought to know by now */
	}
	if (n > 10000) break;	
      }
      esl_abc_GuessAlphabet(ct, ret_type);
    }

  if (*ret_type == eslUNKNOWN) return eslEAMBIGUOUS;
  else                         return eslOK;
}


/* Function:  esl_msa_CreateDigital()
 * Synopsis:  Create a digital <ESL_MSA>.
 * Incept:    SRE, Sun Aug 27 16:49:58 2006 [Leesburg]
 *
 * Purpose:   Same as <esl_msa_Create()>, except the returned MSA is configured
 *            for a digital alignment using internal alphabet <abc>, instead of 
 *            a text alignment.
 *   
 *            Internally, this means the <ax> field is allocated instead of
 *            the <aseq> field, and the <eslMSA_DIGITAL> flag is raised.
 *
 * Args:     <nseq> - number of sequences, or nseq allocation blocksize
 *           <alen> - length of alignment in columns, or -1
 *
 * Returns:   pointer to new MSA object, w/ all values initialized.
 *            Note that <msa->nseq> is initialized to 0, even though space
 *            is allocated.
 *           
 * Throws:    NULL on allocation failure.          
 *
 * Xref:      squid's MSAAlloc()
 */
ESL_MSA *
esl_msa_CreateDigital(const ESL_ALPHABET *abc, int nseq, int64_t alen)
{
  int      status;
  ESL_MSA *msa;
  int      i;

  msa = create_mostly(nseq, alen); /* aseq is null upon successful return */
  if (msa == NULL) return NULL; /* already threw error in mostly_create, so percolate */

  ESL_ALLOC(msa->ax,   sizeof(ESL_DSQ *) * msa->sqalloc); 
  for (i = 0; i < msa->sqalloc; i++)
    msa->ax[i] = NULL;

  if (alen != -1)
    {
      for (i = 0; i < nseq; i++) {
	ESL_ALLOC(msa->ax[i], sizeof(ESL_DSQ) * (alen+2));
	msa->ax[i][0] = msa->ax[i][alen+1] = eslDSQ_SENTINEL; /* help the poor */
      }
      msa->nseq = nseq;
    }

  msa->abc    = (ESL_ALPHABET *) abc; /* this cast away from const-ness is deliberate & safe. */
  msa->flags |= eslMSA_DIGITAL;
  return msa;

 ERROR:
  esl_msa_Destroy(msa);
  return NULL;
}

/* Function:  esl_msa_Digitize()
 * Synopsis:  Digitizes an msa, converting it from text mode.
 * Incept:    SRE, Sat Aug 26 17:33:08 2006 [AA 5302 to Dulles]
 *
 * Purpose:   Given an alignment <msa> in text mode, convert it to
 *            digital mode, using alphabet <abc>.
 *            
 *            Internally, the <ax> digital alignment field is filled,
 *            the <aseq> text alignment field is destroyed and free'd,
 *            a copy of the alphabet pointer is kept in the msa's
 *            <abc> reference, and the <eslMSA_DIGITAL> flag is raised
 *            in <flags>.
 *
 * Args:      abc    - digital alphabet
 *            msa    - multiple alignment to digitize
 *
 * Returns:   <eslOK> on success;
 *            <eslEINVAL> if one or more sequences contain invalid characters
 *            that can't be digitized. If this happens, the <msa> is returned
 *            unaltered - left in text mode, with <aseq> as it was. (This is
 *            a normal error, because <msa->aseq> may be user input that we 
 *            haven't validated yet.)
 *
 * Throws:    <eslEMEM> on allocation failure; in this case, state of <msa> may be 
 *            wedged, and it should only be destroyed, not used.
 */
int
esl_msa_Digitize(const ESL_ALPHABET *abc, ESL_MSA *msa)
{
  int status;
  int i;

  /* Contract checks
   */
  if (msa->aseq == NULL)           ESL_EXCEPTION(eslEINVAL, "msa has no text alignment");
  if (msa->ax   != NULL)           ESL_EXCEPTION(eslEINVAL, "msa already has digital alignment");
  if (msa->flags & eslMSA_DIGITAL) ESL_EXCEPTION(eslEINVAL, "msa is flagged as digital");

  /* Validate before we convert. Then we can leave the <aseq> untouched if
   * any of the sequences contain invalid characters.
   */
  for (i = 0; i < msa->nseq; i++)
    if (esl_abc_ValidateSeq(abc, msa->aseq[i], msa->alen, NULL) != eslOK) 
      return eslEINVAL;

  /* Convert, sequence-by-sequence, free'ing aseq as we go.
   */
  ESL_ALLOC(msa->ax, msa->sqalloc * sizeof(ESL_DSQ *));
  for (i = 0; i < msa->nseq; i++)
    {
      ESL_ALLOC(msa->ax[i], (msa->alen+2) * sizeof(ESL_DSQ));
      status = esl_abc_Digitize(abc, msa->aseq[i], msa->ax[i]);
      if (status != eslOK) goto ERROR;
      free(msa->aseq[i]);
    }    
  for (; i < msa->sqalloc; i++) 
    msa->ax[i] = NULL;
  free(msa->aseq);
  msa->aseq = NULL;

  msa->abc   =  (ESL_ALPHABET *) abc; /* convince compiler that removing const-ness is safe */
  msa->flags |= eslMSA_DIGITAL;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_msa_Textize()
 * Synopsis:  Convert a digital msa to text mode.
 * Incept:    SRE, Sat Aug 26 18:14:30 2006 [AA 5302 to Dulles]
 *
 * Purpose:   Given an alignment <msa> in digital mode, convert it
 *            to text mode.
 *            
 *            Internally, the <aseq> text alignment field is filled, the
 *            <ax> digital alignment field is destroyed and free'd, the
 *            msa's <abc> digital alphabet reference is nullified, and 
 *            the <eslMSA_DIGITAL> flag is dropped in <flags>.
 *            
 * Args:      msa   - multiple alignment to convert to text
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslECORRUPT> if one or more of the digitized alignment strings
 *            contain invalid characters.
 */
int
esl_msa_Textize(ESL_MSA *msa)
{
  int status;
  int i;

  /* Contract checks
   */
  if (msa->ax   == NULL)               ESL_EXCEPTION(eslEINVAL, "msa has no digital alignment");
  if (msa->aseq != NULL)               ESL_EXCEPTION(eslEINVAL, "msa already has text alignment");
  if (! (msa->flags & eslMSA_DIGITAL)) ESL_EXCEPTION(eslEINVAL, "msa is not flagged as digital");
  if (msa->abc  == NULL)               ESL_EXCEPTION(eslEINVAL, "msa has no digital alphabet");

  /* Convert, sequence-by-sequence, free'ing ax as we go.
   */
  ESL_ALLOC(msa->aseq, msa->sqalloc * sizeof(char *));
  for (i = 0; i < msa->nseq; i++)
    {
      ESL_ALLOC(msa->aseq[i], (msa->alen+1) * sizeof(char));
      status = esl_abc_Textize(msa->abc, msa->ax[i], msa->alen, msa->aseq[i]);
      if (status != eslOK) goto ERROR;
      free(msa->ax[i]);
    }
  for (; i < msa->sqalloc; i++)
    msa->aseq[i] = NULL;
  free(msa->ax);
  msa->ax = NULL;
  
  msa->abc    = NULL;      	 /* nullify reference (caller still owns real abc) */
  msa->flags &= ~eslMSA_DIGITAL; /* drop the flag */
  return eslOK;

 ERROR:
  return status;
}



/* Function:  esl_msafile_GuessAlphabet()
 * Synopsis:  Guess what kind of sequences the alignment file contains.
 * Incept:    SRE, Wed May 16 10:45:37 2007 [Janelia]
 *
 * Purpose:   Guess the alphabet of the sequences in the open
 *            <ESL_MSAFILE> -- <eslDNA>, <eslRNA>, or <eslAMINO> --
 *            based on the composition of the next MSA in the
 *            file. Usually this would be the first MSA, because we
 *            would call <esl_msafile_GuessAlphabet()> immediately
 *            after opening a new MSA file.
 *            
 * Returns:   Returns <eslOK> on success, and <*ret_type> is set
 *            to <eslDNA>, <eslRNA>, or <eslAMINO>. 
 *            
 *            Returns <eslEAMBIGUOUS> and sets <*ret_type> to
 *            <eslUNKNOWN> if the first alignment in the file contains
 *            no more than ten residues total, or if its alphabet
 *            cannot be reliably guessed (it contains IUPAC degeneracy
 *            codes, but no amino acid specific residues).
 * 
 *            Returns <eslEFORMAT> if a parse error is encountered
 *            in trying to read the alignment file. <msafp->errbuf>
 *            is set to a useful error message if this occurs; 
 *            <*ret_type> is set to <eslUNKNOWN>.
 *
 *            Returns <eslENODATA> if the file is empty and no
 *            alignment was found; <*ret_type> is set to <eslUNKNOWN>.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Xref:      J1/62.
 */
int
esl_msafile_GuessAlphabet(ESL_MSAFILE *msafp, int *ret_type)
{
  int      status;
  ESL_MSA *msa = NULL;
  
  /* If the msafp is digital mode, we already know the type (so why
   * are we being called?) If do_digital is TRUE, msafp->abc is
   * non-NULL: 
   */
  if (msafp->abc != NULL) { *ret_type = msafp->abc->type; return eslOK; } /* that was easy */

  /* If there's already an MSA cached, we've already called
   * GuessAlphabet(); don't read another one, or we'll overwrite the
   * first.
   */
  if (msafp->msa_cache != NULL) return esl_msa_GuessAlphabet(msafp->msa_cache, ret_type);

  /* Read first alignment, collect its residue composition for
   * passing off to esl_abc_GuessAlphabet()
   */
  status = esl_msa_Read(msafp, &msa);
  if      (status == eslEOF)     return eslENODATA;
  else if (status != eslOK)      return status;

  /* Store the msa in the MSAFILE's cache.  This is because if we're
   * reading from stdin or a gzip pipe, we'll have trouble rewinding
   * to prepare for the first msa_Read() call.
   */
  msafp->msa_cache = msa;

  /* And over to msa_GuessAlphabet() for the decision.
   */
  return esl_msa_GuessAlphabet(msa, ret_type);
}



/* Function:  esl_msafile_OpenDigital()
 * Synopsis:  Open an msa file for digital input.
 * Incept:    SRE, Sun Aug 27 17:40:33 2006 [Leesburg]
 *
 * Purpose:   Same as <esl_msafile_Open()>, except the alignment file
 *            will be read into a digitized internal representation,
 *            using internal alphabet <abc>, rather than the default
 *            internal ASCII text representation.
 *            
 * Args:      abc      - pointer to internal alphabet
 *            filename - name of alignment data file to open;
 *                       if "*.gz", attempt to read through <gunzip -dc> using <popen()>;
 *                       or "-" for stdin 
 *            format   - file format code (e.g. <eslMSAFILE_STOCKHOLM>);
 *                       or <eslMSAFILE_UNKNOWN> to invoke format autodetection.
 *            env      - NULL, or the name of an environment variable from which
 *                       to retrieve a colon-delimited directory list to search
 *                       for <filename> in. (e.g. "HMMERDB")
 *            ret_msafp - RETURN: open MSAFILE.
 *
 * Returns:  <eslOK> on success, and <ret_msafp> is set to point at
 *           an open <ESL_MSAFILE>. Caller frees this file pointer with
 *           <esl_msafile_Close()>.
 *           
 *           <eslENOTFOUND> if <filename> cannot be opened;
 *           <eslEFORMAT> if autodetection is attempted and format
 *           cannot be determined.
 *           
 * Throws:   <eslEMEM> on allocation failure.
 *           <eslEINVAL> if format autodetection is attempted on 
 *           stdin or a gunzip pipe.
 */
int
esl_msafile_OpenDigital(const ESL_ALPHABET *abc, const char *filename, 
			int format, const char *env, ESL_MSAFILE **ret_msafp)
{
  ESL_MSAFILE *msafp;
  int          status;

  if ((status = msafile_open(filename, format, env, &msafp)) != eslOK) return status;
  esl_msafile_SetDigital(msafp, abc);

  *ret_msafp = msafp;
  return eslOK;
}



/* Function:  esl_msafile_SetDigital()
 * Synopsis:  Set an open <ESL_MSAFILE> to read in digital mode.
 * Incept:    SRE, Wed May 16 10:40:24 2007 [Janelia]
 *
 * Purpose:   Given an open <ESL_MSAFILE>, set it so that all subsequent
 *            calls to <esl_msa_Read()> will read multiple alignments
 *            in digital mode instead of text mode, using alphabet
 *            <abc>.
 *            
 *            <esl_msafile_Open(); esl_msafile_SetDigital()> is
 *            equivalent to <esl_msafile_OpenDigital()>. The two-step
 *            version is useful if you don't already know the alphabet
 *            type for your msa file, and you need to call
 *            <esl_msafile_GuessAlphabet()> after opening the file but
 *            before setting its digital alphabet.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_msafile_SetDigital(ESL_MSAFILE *msafp, const ESL_ALPHABET *abc)
{
  msafp->abc        = abc;
  msafp->do_digital = TRUE;
  return eslOK;
}
#endif /* eslAUGMENT_ALPHABET */
/*---------------------- end of digital MSA functions -----------------------*/





/******************************************************************************
 *# 4. Random MSA database access (augmentation: ssi)
 *****************************************************************************/
#ifdef eslAUGMENT_SSI
/* Function:  esl_msafile_PositionByKey()
 * Synopsis:  Use SSI to reposition file to start of named MSA.
 * Incept:    SRE, Mon May 28 11:04:59 2007 [Janelia]
 *
 * Purpose:   Reposition <afp> so that the next MSA we read
 *            will be the one named (or accessioned) <key>.
 *
 * Returns:   <eslOK> on success, and the file <afp> is repositioned
 *            such that the next <esl_msafile_Read()> call will read the
 *            alignment named <key>.
 *            
 *            Returns <eslENOTFOUND> if <key> isn't found in the index
 *            for <afp>. 
 *            
 *            Returns <eslEFORMAT> if something goes wrong trying to
 *            read the index, indicating some sort of file format
 *            problem in the SSI file.
 *
 * Throws:    <eslEMEM> on allocation failure;
 *            <eslEINVAL> if there's no open SSI index;
 *            <eslESYS> if an <fseek()> fails.
 *            In all these cases, the state of the <afp> is uncertain
 *            and may be corrupt; the application should not continue
 *            to use it.
 */
int
esl_msafile_PositionByKey(ESL_MSAFILE *afp, const char *key)
{
  uint16_t fh;
  off_t    offset;
  int      status;

  if (afp->ssi == NULL) ESL_EXCEPTION(eslEINVAL, "Need an open SSI index to call esl_msafile_PositionByKey()");
  if ((status = esl_ssi_FindName(afp->ssi, key, &fh, &offset, NULL, NULL)) != eslOK) return status;
  if (fseeko(afp->f, offset, SEEK_SET) != 0)    ESL_EXCEPTION(eslESYS, "fseek failed");
  
  /* If the <afp> had an MSA cached, we will probably have to discard
   * it, unless by chance it's exactly the MSA we're looking for.
   */
  if (afp->msa_cache != NULL)
    {
      if ( (afp->msa_cache->name == NULL || strcmp(afp->msa_cache->name, key) != 0) &&
	   (afp->msa_cache->acc  == NULL || strcmp(afp->msa_cache->acc,  key) != 0))
	{
	  esl_msa_Destroy(afp->msa_cache);
	  afp->msa_cache = NULL;
	}
    }

  /* The linenumber gets messed up after a file positioning. Best we can do
   * is to reset it to zero.
   */
  afp->linenumber = 0; 
  return eslOK;
}
#endif /*eslAUGMENT_SSI*/
/*------------- end of functions added by SSI augmentation -------------------*/




/******************************************************************************
 *# 5. General i/o API, all alignment formats                                 
 *****************************************************************************/

/* msafile_getline():
 * load the next line of <afp> into <afp->buf>. 
 * Returns eslOK on success, eslEOF on normal eof.
 * Throws eslEMEM on alloc failure.
 */
static int
msafile_getline(ESL_MSAFILE *afp)
{
  int status;
  status = esl_fgets(&(afp->buf), &(afp->buflen), afp->f);
  afp->linenumber++;
  return status;
}

/* maxwidth()
 * Return the length of the longest string in 
 * an array of strings.
 */
static int64_t
maxwidth(char **s, int n)
{
  int64_t max,len;
  int     i; 
  
  max = 0;
  for (i = 0; i < n; i++)
    if (s[i] != NULL)
      {
	len = strlen(s[i]);
	if (len > max) max = len;
      }
  return max;
}

/* actually_write_stockholm()
 * SRE, Fri May 21 17:39:22 1999 [St. Louis]
 *
 * Write an alignment in Stockholm format to an open file. This is the
 * function that actually does the work. The API's WriteStockholm()
 * and WriteStockholmOneBlock() are wrappers.
 *
 * Args:     fp    - file that's open for writing
 *           msa   - alignment to write        
 *           cpl   - characters to write per line in alignment block
 *
 * Returns:  eslOK on success.
 * 
 * Throws:   eslEMEM on allocation failure.
 */
static int
actually_write_stockholm(FILE *fp, const ESL_MSA *msa, int cpl)
{
  int  i, j;
  int  maxname;		/* maximum name length     */
  int  maxgf;		/* max #=GF tag length     */
  int  maxgc;		/* max #=GC tag length     */
  int  maxgr; 		/* max #=GR tag length     */
  int  margin;        	/* total left margin width */
  int  gslen;		/* length of a #=GS tag    */
  char *buf = NULL;
  int  currpos;
  char *s, *tok;
  int  acpl;            /* actual number of character per line */
  int  status;
  
  /* Figure out how much space we need for name + markup
   * to keep the alignment in register. Required by Stockholm
   * spec, even though our Stockholm parser doesn't care (Erik's does).
   *
   * The left margin of an alignment block can be composed of:
   * 
   * <seqname>                      max length: maxname + 1
   * #=GC <gc_tag>                  max length: 4 + 1 + maxgc + 1
   * #=GR <seqname> <gr_tag>        max length: 4 + 1 + maxname + 1 + maxgr + 1
   * 
   * <margin> is the max of these. It is the total length of the
   * left margin that we need to leave, inclusive of the last space.
   * 
   * Then when we output, we do:
   * name:  <leftmargin-1>
   * gc:    #=GC <leftmargin-6>
   * gr:    #=GR <maxname> <leftmargin-maxname-7>
   *
   * xref STL9/p17
   */
  maxname = maxwidth(msa->sqname, msa->nseq);
  
  maxgf   = maxwidth(msa->gf_tag, msa->ngf);
  if (maxgf < 2) maxgf = 2;

  maxgc   = maxwidth(msa->gc_tag, msa->ngc);
  if (msa->rf      !=NULL && maxgc < 2) maxgc = 2;
  if (msa->ss_cons !=NULL && maxgc < 7) maxgc = 7;
  if (msa->sa_cons !=NULL && maxgc < 7) maxgc = 7;

  maxgr   = maxwidth(msa->gr_tag, msa->ngr);
  if (msa->ss != NULL && maxgr < 2) maxgr = 2;
  if (msa->sa != NULL && maxgr < 2) maxgr = 2;

  margin = maxname + 1;
  if (maxgc > 0 && maxgc+6 > margin) margin = maxgc+6;
  if (maxgr > 0 && maxname+maxgr+7 > margin) margin = maxname+maxgr+7; 
  
  /* Allocate a tmp buffer to hold sequence chunks in
   */
  ESL_ALLOC(buf, sizeof(char) * (cpl+1));

  /* Magic Stockholm header
   */
  fprintf(fp, "# STOCKHOLM 1.0\n");

  /* Free text comments
   */
  for (i = 0;  i < msa->ncomment; i++)
    fprintf(fp, "# %s\n", msa->comment[i]);
  if (msa->ncomment > 0) fprintf(fp, "\n");

  /* GF section: per-file annotation
   */
  if (msa->name != NULL) fprintf(fp, "#=GF %-*s %s\n", maxgf, "ID", msa->name);
  if (msa->acc  != NULL) fprintf(fp, "#=GF %-*s %s\n", maxgf, "AC", msa->acc);
  if (msa->desc != NULL) fprintf(fp, "#=GF %-*s %s\n", maxgf, "DE", msa->desc);
  if (msa->au   != NULL) fprintf(fp, "#=GF %-*s %s\n", maxgf, "AU", msa->au);
  
  /* Thresholds are hacky. Pfam has two. Rfam has one.
   */
  if      (msa->cutset[eslMSA_GA1] && msa->cutset[eslMSA_GA2])
    fprintf(fp, "#=GF %-*s %.1f %.1f\n", 
	    maxgf, "GA", msa->cutoff[eslMSA_GA1], msa->cutoff[eslMSA_GA2]);
  else if (msa->cutset[eslMSA_GA1])
    fprintf(fp, "#=GF %-*s %.1f\n", 
	    maxgf, "GA", msa->cutoff[eslMSA_GA1]);

  if      (msa->cutset[eslMSA_NC1] && msa->cutset[eslMSA_NC2])
    fprintf(fp, "#=GF %-*s %.1f %.1f\n",
	    maxgf, "NC", msa->cutoff[eslMSA_NC1], msa->cutoff[eslMSA_NC2]);
  else if (msa->cutset[eslMSA_NC1])
    fprintf(fp, "#=GF %-*s %.1f\n",
	    maxgf, "NC", msa->cutoff[eslMSA_NC1]);

  if      (msa->cutset[eslMSA_TC1] && msa->cutset[eslMSA_TC2])
    fprintf(fp, "#=GF %-*s %.1f %.1f\n",
	    maxgf, "TC", msa->cutoff[eslMSA_TC1], msa->cutoff[eslMSA_TC2]);
  else if (msa->cutset[eslMSA_TC1])
    fprintf(fp, "#=GF %-*s %.1f\n", 
	    maxgf, "TC", msa->cutoff[eslMSA_TC1]);

  for (i = 0; i < msa->ngf; i++)
    fprintf(fp, "#=GF %-*s %s\n", maxgf, msa->gf_tag[i], msa->gf[i]); 
  fprintf(fp, "\n");


  /* GS section: per-sequence annotation
   */
  if (msa->flags & eslMSA_HASWGTS) 
    {
      for (i = 0; i < msa->nseq; i++) 
	fprintf(fp, "#=GS %-*s WT %.2f\n", 
		maxname, msa->sqname[i], msa->wgt[i]);		
      fprintf(fp, "\n");
    }

  if (msa->sqacc != NULL) 
    {
      for (i = 0; i < msa->nseq; i++) 
	if (msa->sqacc[i] != NULL)
	  fprintf(fp, "#=GS %-*s AC %s\n", 
		  maxname, msa->sqname[i], msa->sqacc[i]);
      fprintf(fp, "\n");
    }

  if (msa->sqdesc != NULL) 
    {
      for (i = 0; i < msa->nseq; i++) 
	if (msa->sqdesc[i] != NULL)
	  fprintf(fp, "#=GS %-*s DE %s\n", 
		  maxname, msa->sqname[i], msa->sqdesc[i]);
      fprintf(fp, "\n");
    }

  for (i = 0; i < msa->ngs; i++)
    {
      /* Multiannotated GS tags are possible; for example, 
       *     #=GS foo DR PDB; 1xxx;
       *     #=GS foo DR PDB; 2yyy;
       * These are stored, for example, as:
       *     msa->gs[0][0] = "PDB; 1xxx;\nPDB; 2yyy;"
       * and must be decomposed.
       */
      gslen = strlen(msa->gs_tag[i]);
      for (j = 0; j < msa->nseq; j++)
	if (msa->gs[i][j] != NULL)
	  {
	    s = msa->gs[i][j];
	    while (esl_strtok(&s, "\n", &tok, NULL) == eslOK)
	      fprintf(fp, "#=GS %-*s %-*s %s\n", 
		      maxname, msa->sqname[j],
		      gslen,   msa->gs_tag[i], 
		      tok);
	  }
      fprintf(fp, "\n");
    }

  /* Alignment section:
   * contains aligned sequence, #=GR annotation, and #=GC annotation
   */
  for (currpos = 0; currpos < msa->alen; currpos += cpl)
    {
      acpl = (msa->alen - currpos > cpl)? cpl : msa->alen - currpos;

      if (currpos > 0) fprintf(fp, "\n");
      for (i = 0; i < msa->nseq; i++)
	{
#ifdef eslAUGMENT_ALPHABET
	  if (msa->flags & eslMSA_DIGITAL)
	    esl_abc_TextizeN(msa->abc, msa->ax[i] + currpos + 1, acpl, buf);
	  else
	    strncpy(buf, msa->aseq[i] + currpos, acpl);
#else
	  strncpy(buf, msa->aseq[i] + currpos, acpl);
#endif
	  
	  buf[acpl] = '\0';	      
	  fprintf(fp, "%-*s %s\n", 
		  margin-1, msa->sqname[i], buf);

	  if (msa->ss != NULL && msa->ss[i] != NULL) {
	    strncpy(buf, msa->ss[i] + currpos, acpl);
	    buf[acpl] = '\0';	 
	    fprintf(fp, "#=GR %-*s %-*s %s\n", 
		    maxname,          msa->sqname[i],
		    margin-maxname-7, "SS",
		    buf);
	  }
	  if (msa->sa != NULL && msa->sa[i] != NULL) {
	    strncpy(buf, msa->sa[i] + currpos, acpl);
	    buf[acpl] = '\0';
	    fprintf(fp, "#=GR %-*s %-*s %s\n",
		    maxname,          msa->sqname[i],
		    margin-maxname-7, "SA",
		    buf);
	  }
	  for (j = 0; j < msa->ngr; j++)
	    if (msa->gr[j][i] != NULL) {
	      strncpy(buf, msa->gr[j][i] + currpos, acpl);
	      buf[acpl] = '\0';
	      fprintf(fp, "#=GR %-*s %-*s %s\n", 
		      maxname,          msa->sqname[i],
		      margin-maxname-7, msa->gr_tag[j],
		      buf);
	    }
	}
      if (msa->ss_cons != NULL) {
	strncpy(buf, msa->ss_cons + currpos, acpl);
	buf[acpl] = '\0';
	fprintf(fp, "#=GC %-*s %s\n", margin-6, "SS_cons", buf);
      }

      if (msa->sa_cons != NULL) {
	strncpy(buf, msa->sa_cons + currpos, acpl);
	buf[acpl] = '\0';
	fprintf(fp, "#=GC %-*s %s\n", margin-6, "SA_cons", buf);
      }

      if (msa->rf != NULL) {
	strncpy(buf, msa->rf + currpos, acpl);
	buf[acpl] = '\0';
	fprintf(fp, "#=GC %-*s %s\n", margin-6, "RF", buf);
      }
      for (j = 0; j < msa->ngc; j++) {
	strncpy(buf, msa->gc[j] + currpos, acpl);
	buf[acpl] = '\0';
	fprintf(fp, "#=GC %-*s %s\n", margin-6, msa->gc_tag[j], buf);
      }
    }
  fprintf(fp, "//\n");
  free(buf);
  return eslOK;

 ERROR:
  if (buf != NULL) free(buf);
  return status;
}


/* write_stockholm():
 * SRE, Fri Jan 28 09:24:02 2005 [St. Louis]
 *
 * Purpose:   Write an alignment <msa> in Stockholm format 
 *            to a stream <fp>, in multiblock format, with
 *            50 residues per line.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      squid's WriteStockholm(), 1999.
 */
static int
write_stockholm(FILE *fp, const ESL_MSA *msa)
{
  return (actually_write_stockholm(fp, msa, 50)); /* 50 char per block */
}

/* write_pfam():
 * SRE, Fri Jan 28 09:25:42 2005 [St. Louis]
 *
 * Purpose:   Write an alignment <msa> in Stockholm format 
 *            to a stream <fp>, in single block (Pfam) format.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      squid's WriteStockholmOneBlock(), 1999.
 */
static int
write_pfam(FILE *fp, const ESL_MSA *msa)
{
  return (actually_write_stockholm(fp, msa, msa->alen)); /* one big block */
}


static int
is_blankline(char *s)
{
  for (; *s != '\0'; s++)
    if (! isspace((int) *s)) return FALSE;
  return TRUE;
}

/* Format of a GF line:
 *    #=GF <tag> <text>
 * Returns eslOK on success; eslEFORMAT on parse failure.
 * Throws eslEMEM on allocation failure.
 */
static int
parse_gf(ESL_MSA *msa, char *buf)
{
  char *gf;
  char *tag;
  char *text;
  char *tok;
  char *s;
  int   n;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &gf,   NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag,  NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, "\n\r",    &text, &n)   != eslOK) return eslEFORMAT;
  while (*text && (*text == ' ' || *text == '\t')) text++;

  if      (strcmp(tag, "ID") == 0)
    status = esl_strdup(text, n, &(msa->name));
  else if (strcmp(tag, "AC") == 0) 
    status = esl_strdup(text, n, &(msa->acc));
  else if (strcmp(tag, "DE") == 0) 
    status = esl_strdup(text, n, &(msa->desc));
  else if (strcmp(tag, "AU") == 0) 
    status = esl_strdup(text, n, &(msa->au));
  else if (strcmp(tag, "GA") == 0) 
    {				/* Pfam has GA1, GA2. Rfam just has GA1. */
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &tok, NULL)) != eslOK) 
	return eslEFORMAT;
      msa->cutoff[eslMSA_GA1] = atof(tok);
      msa->cutset[eslMSA_GA1] = TRUE;
      if ((esl_strtok(&s, " \t\n\r", &tok, NULL)) == eslOK) 
	{
	  msa->cutoff[eslMSA_GA2] = atof(tok);
	  msa->cutset[eslMSA_GA2] = TRUE;
	}
      status = eslOK;
    }
  else if (strcmp(tag, "NC") == 0) 
    {
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &tok, NULL)) != eslOK) 
	return eslEFORMAT;
      msa->cutoff[eslMSA_NC1] = atof(tok);
      msa->cutset[eslMSA_NC1] = TRUE;
      if ((esl_strtok(&s, " \t\n\r", &tok, NULL)) == eslOK) 
	{
	  msa->cutoff[eslMSA_NC2] = atof(tok);
	  msa->cutset[eslMSA_NC2] = TRUE;
	}
      status = eslOK;
    }
  else if (strcmp(tag, "TC") == 0) 
    {
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &tok, NULL)) != eslOK) 
	return eslEFORMAT;
      msa->cutoff[eslMSA_TC1] = atof(tok);
      msa->cutset[eslMSA_TC1] = TRUE;
      if ((esl_strtok(&s, "\t\n\r", &tok, NULL)) == eslOK) 
	{
	  msa->cutoff[eslMSA_TC2] = atof(tok);
	  msa->cutset[eslMSA_TC2] = TRUE;
	}
      status = eslOK;
    }
  else 				/* an unparsed #=GF: */
    status = esl_msa_AddGF(msa, tag, text);

  return status;
}


/* Format of a GS line:
 *    #=GS <seqname> <tag> <text>
 * Return <eslOK> on success; <eslEFORMAT> on parse error.
 * Throws <eslEMEM> on allocation error (trying to grow for a new
 *        name; <eslEINVAL> if we try to grow an ungrowable MSA.
 */
static int
parse_gs(ESL_MSA *msa, char *buf)
{
  char *gs;
  char *seqname;
  char *tag;
  char *text; 
  int   seqidx;
  char *s;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &gs,      NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &seqname, NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag,     NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, "\n\r",    &text,    NULL) != eslOK) return eslEFORMAT;
  while (*text && (*text == ' ' || *text == '\t')) text++;
  
  /* GS usually follows another GS; guess lastidx+1 */
  status = get_seqidx(msa, seqname, msa->lastidx+1, &seqidx);
  if (status != eslOK) return status;
  msa->lastidx = seqidx;

  if (strcmp(tag, "WT") == 0)
    {
      msa->wgt[seqidx] = atof(text);
      msa->flags      |= eslMSA_HASWGTS;
      status           = eslOK;
    }
  else if (strcmp(tag, "AC") == 0)
    status = set_seq_accession(msa, seqidx, text);
  else if (strcmp(tag, "DE") == 0)
    status = set_seq_description(msa, seqidx, text);
  else				
    status = esl_msa_AddGS(msa, tag, seqidx, text);

  return status;
}


/* parse_gc():
 * Format of a GC line:
 *    #=GC <tag> <aligned text>
 */
static int 
parse_gc(ESL_MSA *msa, char *buf)
{
  char *gc;
  char *tag;
  char *text; 
  char *s;
  int   len;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &gc,   NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag,  NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &text, &len) != eslOK) return eslEFORMAT;
  
  if (strcmp(tag, "SS_cons") == 0)
    status = esl_strcat(&(msa->ss_cons), -1, text, len);
  else if (strcmp(tag, "SA_cons") == 0)
    status = esl_strcat(&(msa->sa_cons), -1, text, len);
  else if (strcmp(tag, "RF") == 0)
    status = esl_strcat(&(msa->rf), -1, text, len);
  else
    status = esl_msa_AppendGC(msa, tag, text);

  return status;
}

/* parse_gr():
 * Format of a GR line:
 *    #=GR <seqname> <featurename> <text>
 */
static int
parse_gr(ESL_MSA *msa, char *buf)
{
  char *gr;
  char *seqname;
  char *tag;
  char *text;
  int   seqidx;
  int   len;
  int   j;
  char *s;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &gr,      NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &seqname, NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag,     NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &text,    &len) != eslOK) return eslEFORMAT;

  /* GR usually follows sequence it refers to; guess msa->lastidx */
  status = get_seqidx(msa, seqname, msa->lastidx, &seqidx);
  if (status != eslOK) return status;
  msa->lastidx = seqidx;

  if (strcmp(tag, "SS") == 0) 
    {
      if (msa->ss == NULL)
	{
	  ESL_ALLOC(msa->ss,    sizeof(char *) * msa->sqalloc);
	  ESL_ALLOC(msa->sslen, sizeof(int64_t)* msa->sqalloc);
	  for (j = 0; j < msa->sqalloc; j++)
	    {
	      msa->ss[j]    = NULL;
	      msa->sslen[j] = 0;
	    }
	}
      status = esl_strcat(&(msa->ss[seqidx]), msa->sslen[seqidx], text, len);
      msa->sslen[seqidx] += len;
    }
  else if (strcmp(tag, "SA") == 0)
    {
      if (msa->sa == NULL)
	{
	  ESL_ALLOC(msa->sa,    sizeof(char *) * msa->sqalloc);
	  ESL_ALLOC(msa->salen, sizeof(int64_t)* msa->sqalloc);
	  for (j = 0; j < msa->sqalloc; j++) 
	    {
	      msa->sa[j]    = NULL;
	      msa->salen[j] = 0;
	    }
	}
      status = esl_strcat(&(msa->sa[seqidx]), msa->salen[seqidx], text, len);
      msa->salen[seqidx] += len;
    }
  else 
    status = esl_msa_AppendGR(msa, tag, seqidx, text);
  return status;

 ERROR:
  return status;
}


/* parse_comment():
 * comments are simply stored verbatim, not parsed
 */
static int
parse_comment(ESL_MSA *msa, char *buf)
{
  char *s;
  char *comment;

  s = buf + 1;			               /* skip leading '#' */
  if (*s == '\n' || *s == '\r') { *s = '\0'; comment = s; }  /* deal with blank comment */
  else if (esl_strtok(&s, "\n\r", &comment, NULL)!= eslOK) return eslEFORMAT;
  return (esl_msa_AddComment(msa, comment));
}

/* parse_sequence():
 * Format of line is:
 *     <name>  <aligned text>
 * 
 * On digital sequence, returns <eslEINVAL> if any of the residues can't be digitized.
 */
static int
parse_sequence(ESL_MSA *msa, char *buf)
{
  char *s;
  char *seqname;
  char *text;
  int   seqidx;
  int   len;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &seqname, NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &text,    &len) != eslOK) return eslEFORMAT; 
  
  /* seq usually follows another seq; guess msa->lastidx +1 */
  status = get_seqidx(msa, seqname, msa->lastidx+1, &seqidx);
  if (status != eslOK) return status;
  msa->lastidx = seqidx;

#ifdef eslAUGMENT_ALPHABET
  if (msa->flags & eslMSA_DIGITAL)
    {
      status = esl_abc_dsqcat(msa->abc, &(msa->ax[seqidx]), &(msa->sqlen[seqidx]), text, len);
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL))
    {
      status = esl_strcat(&(msa->aseq[seqidx]), msa->sqlen[seqidx], text, len);
      msa->sqlen[seqidx] += len;
    }

  return status;
}

/* read_stockholm():
 * SRE, Sun Jan 23 08:33:32 2005 [St. Louis]
 *
 * Purpose:   Parse the next alignment from an open Stockholm format alignment
 *            file <afp>, leaving the alignment in <ret_msa>.
 *
 * Returns:   <eslOK> on success, and the alignment is in <ret_msa>.
 *            Returns <eslEOF> if there are no more alignments in <afp>,
 *            and <ret_msa> is set to NULL.
 *            <eslEFORMAT> if parse fails because of a file format problem,
 *            in which case afp->errbuf is set to contain a formatted message 
 *            that indicates the cause of the problem, and <ret_msa> is
 *            set to NULL. 
 *
 *            Returns <eslEINVAL> if we're trying to read a digital alignment,
 *            and an invalid residue is found that can't be digitized.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Xref:      squid's ReadStockholm(), 1999.
 */
static int
read_stockholm(ESL_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA   *msa = NULL;
  char      *s;
  int        status;
  int        status2;
#ifdef eslAUGMENT_SSI
  off_t      offset;
#endif

  if (feof(afp->f))  { status = eslEOF; goto ERROR; }
  afp->errbuf[0] = '\0';

  /* Initialize allocation of the MSA:
   * make it growable, by giving it an initial blocksize of
   * 16 seqs of 0 length.
   */
#ifdef eslAUGMENT_ALPHABET
  if (afp->do_digital == TRUE &&
      (msa = esl_msa_CreateDigital(afp->abc, 16, -1))  == NULL) 
    { status = eslEMEM; goto ERROR; }

#endif
  if (afp->do_digital == FALSE &&
      (msa = esl_msa_Create(16, -1))  == NULL)
    { status = eslEMEM; goto ERROR; }
  if (msa == NULL)    
    { status = eslEMEM; goto ERROR; }

  /* Check the magic Stockholm header line.
   * We have to skip blank lines here, else we perceive
   * trailing blank lines in a file as a format error when
   * reading in multi-record mode.
   */
  do {
#ifdef eslAUGMENT_SSI
    offset = ftello(afp->f);
#endif
    if ((status = msafile_getline(afp)) != eslOK) goto ERROR;
  } while (is_blankline(afp->buf));

  if (strncmp(afp->buf, "# STOCKHOLM 1.", 14) != 0)
    { 
      sprintf(afp->errbuf, "missing \"# STOCKHOLM\" header");
      status = eslEFORMAT; 
      goto ERROR;
    } 

#ifdef eslAUGMENT_SSI
  msa->offset = offset;
#endif

  /* Read the alignment file one line at a time.
   */
  while ((status2 = msafile_getline(afp)) == eslOK) 
    {
      s = afp->buf;
      while (*s == ' ' || *s == '\t') s++;  /* skip leading whitespace */

      if (*s == '#') {

	if      (strncmp(s, "#=GF", 4) == 0)
	  {
	    if ((status = parse_gf(msa, s)) != eslOK)
	      { sprintf(afp->errbuf, "failed to parse #=GF line"); goto ERROR; }
	  }

	else if (strncmp(s, "#=GS", 4) == 0)
	  {
	    if ((status = parse_gs(msa, s)) != eslOK)
	      {	sprintf(afp->errbuf, "failed to parse #=GS line"); goto ERROR; }
	  }

	else if (strncmp(s, "#=GC", 4) == 0)
	  {
	    if  ((status = parse_gc(msa, s)) != eslOK)
	      {	sprintf(afp->errbuf, "failed to parse #=GC line"); goto ERROR; }
	  }

	else if (strncmp(s, "#=GR", 4) == 0)
	  {
	    if ((status = parse_gr(msa, s)) != eslOK)
	      {	sprintf(afp->errbuf, "failed to parse #=GR line"); goto ERROR; }
	  }

	else if ((status = parse_comment(msa, s)) != eslOK)
	  { sprintf(afp->errbuf, "failed to parse comment line"); goto ERROR; }
      } 
      else if (strncmp(s, "//",   2) == 0)   break; /* normal way out */
      else if (*s == '\n' || *s == '\r')     continue;
      else if ((status = parse_sequence(msa, s)) != eslOK)
	{ sprintf(afp->errbuf, "failed to parse sequence line"); goto ERROR; }
    }
  /* If we saw a normal // end, we would've successfully read a line,
   * so when we get here, status (from the line read) should be eslOK.
   */ 
  if (status2 != eslOK)
    { 
      sprintf(afp->errbuf, "didn't find // at end of alignment %.128s",
	      msa->name == NULL ? "" : msa->name);
      status = eslEFORMAT;
      goto ERROR;
    } 
  
  /* Stockholm fmt is complex, so give the newly parsed MSA a good
   * going-over, and finalize the fields of the MSA data structure.
   * verify_parse will fill in errbuf if it sees a problem.
   */
  if (verify_parse(msa, afp->errbuf) != eslOK)
    { status = eslEFORMAT; goto ERROR; } 

  if (ret_msa != NULL) *ret_msa = msa; else esl_msa_Destroy(msa);
  return eslOK;

 ERROR:
  if (msa != NULL)      esl_msa_Destroy(msa);
  if (ret_msa != NULL) *ret_msa = NULL;
  return status;
}



/* Function:  esl_msa_Read()
 * Synopsis:  Read next MSA from a file.
 * Incept:    SRE, Fri Jan 28 08:10:49 2005 [St. Louis]
 *
 * Purpose:   Reads the next MSA from an open MSA file <afp>,
 *            and returns it via <ret_msa>. 
 *
 * Returns:   <eslOK> on success, and <ret_msa> points at the
 *            new MSA object.
 *            
 *            Returns <eslEOF> if there are no more alignments in the file.
 *            
 *            Returns <eslEFORMAT> if there is a parse error, and <afp->errbuf>
 *            is set to an informative message.
 *            
 *            <eslEINVAL> if we're trying to read a digital alignment,
 *            but one or more residues are seen in the file that
 *            aren't valid in our alphabet.
 *            
 * Throws:    <eslEMEM> on allocation failure.           
 *            <eslEINCONCEIVABLE> on internal error.
 */
int
esl_msa_Read(ESL_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA *msa;
  int      status;

  *ret_msa = NULL;

  /* If we've just used GuessAlphabet(), we have an MSA already read
   * and stored in the MSAFILE's cache. Just return it, after worrying
   * about whether it's supposed to be in digital or text mode. (It
   * should always be in text mode, and maybe in need of Digitize(),
   * given how GuessAlphabet works now; but this is coded for more
   * generality in case we use the MSA cache some other way in the
   * future.)
   */
  if (afp->msa_cache != NULL) 
    {
#ifdef eslAUGMENT_ALPHABET
      if      (afp->do_digital   && !(afp->msa_cache->flags & eslMSA_DIGITAL)) {
	if ((status = esl_msa_Digitize(afp->abc, afp->msa_cache)) != eslOK) return status; 
      }
      else if (! afp->do_digital && (afp->msa_cache->flags & eslMSA_DIGITAL)) {
	if ((status = esl_msa_Textize(afp->msa_cache)) != eslOK) return status;
      }
#endif

      *ret_msa         = afp->msa_cache;
      afp->msa_cache = NULL;
      return eslOK;
    }

  /* Otherwise, read the next MSA from the file.
   */      
  switch (afp->format) {
  case eslMSAFILE_STOCKHOLM: status = read_stockholm(afp, &msa); break;
  case eslMSAFILE_PFAM:      status = read_stockholm(afp, &msa); break;
  default:
    ESL_EXCEPTION(eslEINCONCEIVABLE, "no such format");
  }

  *ret_msa = msa;
  return status;
}

/* Function:  esl_msa_Write()
 * Synopsis:  Write an MSA to a file.
 * Incept:    SRE, Fri Jan 28 09:29:28 2005 [St. Louis]
 *
 * Purpose:   Writes an alignment <msa> to an open stream <fp>,
 *            in format specified by <fmt>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslEINCONCEIVABLE> on internal error.
 */
int
esl_msa_Write(FILE *fp, const ESL_MSA *msa, int fmt)
{
  int status;
  switch (fmt) {
  case eslMSAFILE_STOCKHOLM: status = write_stockholm(fp, msa); break;
  case eslMSAFILE_PFAM:      status = write_pfam(fp, msa);      break;
  default: 
    ESL_EXCEPTION(eslEINCONCEIVABLE, "no such format");
  } 
  return status;
}


/* Function:  esl_msa_DescribeFormat()
 * Synopsis:  Convert internal file format code to text string.
 * Incept:    SRE, Fri May 18 11:59:58 2007 [Janelia]
 *
 * Purpose:   Given an internal file format code <fmt> 
 *            (<eslMSAFILE_STOCKHOLM>, for example), returns
 *            a string suitable for printing ("Stockholm",
 *            for example).
 *            
 * Note:      Keep in sync with <esl_sqio_DescribeFormat()>.
 *            The reason we don't just use <esl_sqio_DescribeFormat()>
 *            is so the msa module can be used without the sqio
 *            module.
 */
char *
esl_msa_DescribeFormat(int fmt)
{
  switch (fmt) {
  case eslMSAFILE_UNKNOWN:   return "unknown";
  case eslMSAFILE_STOCKHOLM: return "Stockholm";
  case eslMSAFILE_PFAM:      return "Pfam";
  default:                   esl_fatal("no such format code");
  }
  /*NOTREACHED*/
  return NULL;
}


/* Function:  esl_msa_GuessFileFormat()
 * Synopsis:  Determine the format of an open MSA file.
 * Incept:    SRE, Fri Jan 28 07:29:00 2005 [St. Louis]
 *
 * Purpose:   Attempts to determine the format of an open alignment file
 *            <afp>, for which <afp->format> is <eslMSAFILE_UNKNOWN>. 
 *            If successful, sets <afp->format>.
 *            
 *            Currently a placeholder: it always guesses Stockholm!
 *
 * Returns:   <eslOK> on success, and sets <afp->format>. 
 *            <eslEFORMAT> if format can't be determined.
 *
 * Xref:      squid's MSAFileFormat()
 */
int
esl_msa_GuessFileFormat(ESL_MSAFILE *afp)
{
  /* Placeholder: FIXME: autodetection code goes here.
   */
  afp->format = eslMSAFILE_STOCKHOLM;
  return eslOK;
}
/*-------------------- end of general i/o functions -------------------------*/






/*****************************************************************
 *# 6. Miscellaneous functions for manipulating MSAs
 *****************************************************************/

/* Function:  esl_msa_SequenceSubset()
 * Synopsis:  Select subset of sequences into a smaller MSA.
 * Incept:    SRE, Wed Apr 13 10:05:44 2005 [St. Louis]
 *
 * Purpose:   Given an array <useme> (0..nseq-1) of TRUE/FALSE flags for each
 *            sequence in an alignment <msa>; create a new alignment containing
 *            only those seqs which are flagged <useme=TRUE>. Return a pointer
 *            to this newly allocated alignment through <ret_new>. Caller is
 *            responsible for freeing it.
 *            
 *            The smaller alignment might now contain columns
 *            consisting entirely of gaps or missing data, depending
 *            on what sequence subset was extracted. The caller may
 *            want to immediately call <esl_msa_MinimGaps()> on the
 *            new alignment to clean this up.
 *
 *            Unparsed GS and GR Stockholm annotation that is presumably still
 *            valid is transferred to the new alignment. Unparsed GC, GF, and
 *            comments that are potentially invalidated by taking the subset
 *            of sequences are not transferred to the new MSA.
 *            
 *            Weights are transferred exactly. If they need to be
 *            renormalized to some new total weight (such as the new,
 *            smaller total sequence number), the caller must do that.
 *            
 *            <msa> may be in text mode or digital mode. The new MSA
 *            in <ret_new> will have the same mode.
 *
 * Returns:   <eslOK> on success, and <ret_new> is set to point at a new
 *            (smaller) alignment.
 *
 * Throws:    <eslEINVAL> if the subset has no sequences in it;
 *            <eslEMEM> on allocation error.
 *
 * Xref:      squid's MSASmallerAlignment(), 1999.
 */
int
esl_msa_SequenceSubset(const ESL_MSA *msa, const int *useme, ESL_MSA **ret_new)
{
  ESL_MSA *new = NULL;
  int  nnew;			/* number of seqs in the new MSA */
  int  oidx, nidx;		/* old, new indices */
  int  i;
  int  status;
  
  *ret_new = NULL;

  nnew = 0; 
  for (oidx = 0; oidx < msa->nseq; oidx++)
    if (useme[oidx]) nnew++;
  if (nnew == 0) ESL_EXCEPTION(eslEINVAL, "No sequences selected");

  /* Note that the Create() calls allocate exact space for the sequences,
   * so we will strcpy()/memcpy() into them below.
   */
#ifdef eslAUGMENT_ALPHABET
  if ((msa->flags & eslMSA_DIGITAL) &&
      (new = esl_msa_CreateDigital(msa->abc, nnew, msa->alen)) == NULL)
    {status = eslEMEM; goto ERROR; }
#endif
  if (! (msa->flags & eslMSA_DIGITAL) &&
      (new = esl_msa_Create(nnew, msa->alen)) == NULL) 
    {status = eslEMEM; goto ERROR; }
  if (new == NULL) 
    {status = eslEMEM; goto ERROR; }
  

  /* Copy the old to the new */
  for (nidx = 0, oidx = 0; oidx < msa->nseq; oidx++)
    if (useme[oidx])
      {
#ifdef eslAUGMENT_ALPHABET
	if (msa->flags & eslMSA_DIGITAL)
	  memcpy(new->ax[nidx], msa->ax[oidx], sizeof(ESL_DSQ) * (msa->alen+2));
#endif
	if (! (msa->flags & eslMSA_DIGITAL))
	  strcpy(new->aseq[nidx], msa->aseq[oidx]);
	if ((status = esl_strdup(msa->sqname[oidx], -1, &(new->sqname[nidx])))    != eslOK) goto ERROR;

	new->wgt[nidx] = msa->wgt[oidx];
      
	if (msa->sqacc != NULL && msa->sqacc[oidx] != NULL) {
	  if ((status = set_seq_accession(new, nidx, msa->sqacc[oidx])) != eslOK) goto ERROR;
	}
	if (msa->sqdesc != NULL && msa->sqdesc[oidx] != NULL) {
	  if ((status = set_seq_description(new, nidx, msa->sqdesc[oidx])) != eslOK) goto ERROR;
	}
	if (msa->ss != NULL && msa->ss[oidx] != NULL) {
	  if ((status = set_seq_ss(new, nidx, msa->ss[oidx])) != eslOK) goto ERROR;
	}
	if (msa->sa != NULL && msa->sa[oidx] != NULL) {
	  if ((status = set_seq_sa(new, nidx, msa->sa[oidx])) != eslOK) goto ERROR;
	}
	/* unparsed annotation */
	for(i = 0; i < msa->ngs; i++) {
	  if(msa->gs[i] != NULL) 
	    if ((status = esl_msa_AddGS(new, msa->gs_tag[i], nidx, msa->gs[i][oidx])) != eslOK) goto ERROR;
	}
	for(i = 0; i < msa->ngr; i++) {
	  if(msa->gr[i] != NULL) 
	    if ((status = esl_msa_AppendGR(new, msa->gr_tag[i], nidx, msa->gr[i][oidx])) != eslOK) goto ERROR;
	}

	nidx++;
      }

  new->flags = msa->flags;

  if ((status = esl_strdup(msa->name, -1, &(new->name))) != eslOK) goto ERROR;
  if ((status = esl_strdup(msa->desc, -1, &(new->desc))) != eslOK) goto ERROR;
  if ((status = esl_strdup(msa->acc,  -1, &(new->acc)))  != eslOK) goto ERROR;
  if ((status = esl_strdup(msa->au,   -1, &(new->au)))   != eslOK) goto ERROR;
  if ((status = esl_strdup(msa->ss_cons, msa->alen, &(new->ss_cons))) != eslOK) goto ERROR;
  if ((status = esl_strdup(msa->sa_cons, msa->alen, &(new->sa_cons))) != eslOK) goto ERROR;
  if ((status = esl_strdup(msa->rf, msa->alen, &(new->rf))) != eslOK) goto ERROR;

  for (i = 0; i < eslMSA_NCUTS; i++) {
    new->cutoff[i] = msa->cutoff[i];
    new->cutset[i] = msa->cutset[i];
  }
  
  new->nseq  = nnew;
  new->sqalloc = nnew;

  /* Since we have a fully constructed MSA, we don't need the
   * aux info used by parsers.
   */
  if (new->sqlen != NULL) { free(new->sqlen);  new->sqlen = NULL; }
  if (new->sslen != NULL) { free(new->sslen);  new->sslen = NULL; }
  if (new->salen != NULL) { free(new->salen);  new->salen = NULL; }
  new->lastidx = -1;

  *ret_new = new;
  return eslOK;

 ERROR:
  if (new != NULL) esl_msa_Destroy(new);
  *ret_new = NULL;
  return status;
}


/* remove_broken_basepairs()
 * 
 * Given an array <useme> (0..alen-1) of TRUE/FALSE flags, remove
 * any basepair from SS_cons and individual SS annotation in alignment
 * columns (i,j) for which either <useme[i-1]> or <useme[j-1]> is FALSE.
 * Called automatically from esl_msa_ColumnSubset() with same <useme>. 
 * 
 * If the original structure data is inconsistent it's left untouched.
 *
 * Returns:   <eslOK> on success.
 *            <eslESYNTAX> if WUSS string for SS_cons or msa->ss 
 *            following esl_wuss_nopseudo() is inconsistent.
 *            <eslEINVAL> if a derived ct array implies a pknotted 
 *            SS, this should be impossible
 */
static int
remove_broken_basepairs(ESL_MSA *msa, char *errbuf, const int *useme)
{
  int status;
  int64_t apos;         /* alignment position */
  int  i;
  int *ct;		/* 0..alen-1 base pair partners array for current sequence */
  char *ss_nopseudo;    /* no-pseudoknot version of structure */
  ESL_ALLOC(ct,  sizeof(int)  * (msa->alen+1));
  ESL_ALLOC(ss_nopseudo, sizeof(char) * (msa->alen+1));

  if (msa->ss_cons != NULL) { 
    esl_wuss_nopseudo(msa->ss_cons, ss_nopseudo);
    if ((status = esl_wuss2ct(ss_nopseudo, msa->alen, ct)) != eslOK) ESL_FAIL(status, errbuf, "Consensus structure string is inconsistent.");
    for (apos = 1; apos <= msa->alen; apos++) { 
      if (!(useme[apos-1])) { 
	if (ct[apos] != 0) ct[ct[apos]] = 0;
	ct[apos] = 0;
      }
      /* convert to WUSS SS string and supplant msa->ss_cons */
      if ((status = esl_ct2wuss(ct, msa->alen, msa->ss_cons)) != eslOK) ESL_FAIL(status, errbuf, "Error converting de-knotted bp ct arry to WUSS notation.");
    }
  }
  /* do the same for per-seq SS annotation */
  if (msa->ss != NULL) { 
    for (i = 0; i < msa->nseq; i++) {
      if (msa->ss[i] != NULL) { 
	esl_wuss_nopseudo(msa->ss[i], ss_nopseudo);
	if ((status = esl_wuss2ct(ss_nopseudo, msa->alen, ct)) != eslOK) ESL_FAIL(status, errbuf, "Secondary structure string for seq %d is inconsistent.", i);
	for (apos = 1; apos <= msa->alen; apos++) { 
	  if (!(useme[apos-1])) { 
	    if (ct[apos] != 0) ct[ct[apos]] = 0;
	    ct[apos] = 0;
	  }
	  /* convert to WUSS SS string and supplant msa->ss[i] */
	  if ((status = esl_ct2wuss(ct, msa->alen, msa->ss[i])) != eslOK) ESL_FAIL(status, errbuf, "Error converting de-knotted bp ct arry to WUSS notation.");
	}
      }
    }
  }
  free(ss_nopseudo);
  free(ct);
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.");
  return status; /* NEVERREACHED */
}  

/* Function:  esl_msa_ColumnSubset()
 * Synopsis:  Remove a selected subset of columns from the MSA
 *
 * Incept:    SRE, Sun Feb 27 10:05:07 2005
 *            From squid's MSAShorterAlignment(), 1999
 * 
 * Purpose:   Given an array <useme> (0..alen-1) of TRUE/FALSE flags,
 *            where TRUE means "keep this column in the new alignment"; 
 *            remove all columns annotated as FALSE in the <useme> 
 *            array. This is done in-place on the MSA, so the MSA is 
 *            modified: <msa->alen> is reduced, <msa->aseq> is shrunk 
 *            (or <msa->ax>, in the case of a digital mode alignment), 
 *            and all associated per-residue or per-column annotation
 *            is shrunk.
 * 
 * Returns:   <eslOK> on success.
 *            Possibilities from remove_broken_basepairs() call:
 *            <eslESYNTAX> if WUSS string for SS_cons or msa->ss 
 *            following esl_wuss_nopseudo() is inconsistent.
 *            <eslEINVAL> if a derived ct array implies a pknotted SS.
 */
int
esl_msa_ColumnSubset(ESL_MSA *msa, char *errbuf, const int *useme)
{
  int     status;
  int64_t opos;			/* position in original alignment */
  int64_t npos;			/* position in new alignment      */
  int     idx;			/* sequence index */
  int     i;			/* markup index */

  /* Remove any basepairs from SS_cons and individual sequence SS
   * for aln columns i,j for which useme[i-1] or useme[j-1] are FALSE 
   */
  if((status = remove_broken_basepairs(msa, errbuf, useme)) != eslOK) return status;

  /* Since we're minimizing, we can overwrite in place, within the msa
   * we've already got. 
   * opos runs all the way to msa->alen to include (and move) the \0
   * string terminators (or sentinel bytes, in the case of digital mode)
   */
  for (opos = 0, npos = 0; opos <= msa->alen; opos++)
    {
      if (opos < msa->alen && useme[opos] == FALSE) continue;

      if (npos != opos)	/* small optimization */
	{
	  /* The alignment, and per-residue annotations */
	  for (idx = 0; idx < msa->nseq; idx++)
	    {
#ifdef eslAUGMENT_ALPHABET
	      if (msa->flags & eslMSA_DIGITAL) /* watch off-by-one in dsq indexing */
		msa->ax[idx][npos+1] = msa->ax[idx][opos+1];
	      else
		msa->aseq[idx][npos] = msa->aseq[idx][opos];
#else
	      msa->aseq[idx][npos] = msa->aseq[idx][opos];
#endif /*eslAUGMENT_ALPHABET*/
	      if (msa->ss != NULL && msa->ss[idx] != NULL)
		msa->ss[idx][npos] = msa->ss[idx][opos];
	      if (msa->sa != NULL && msa->sa[idx] != NULL)
		msa->sa[idx][npos] = msa->sa[idx][opos];
	      for (i = 0; i < msa->ngr; i++)
		if (msa->gr[i][idx] != NULL)
		  msa->gr[i][idx][npos] = msa->gr[i][idx][opos];
	    }	  
	  /* The per-column annotations */
	  if (msa->ss_cons != NULL) msa->ss_cons[npos] = msa->ss_cons[opos];
	  if (msa->sa_cons != NULL) msa->sa_cons[npos] = msa->sa_cons[opos];
	  if (msa->rf      != NULL) msa->rf[npos]      = msa->rf[opos];
	  for (i = 0; i < msa->ngc; i++)
	    msa->gc[i][npos] = msa->gc[i][opos];
	}
      npos++;
    }
  msa->alen = npos-1;	/* -1 because npos includes NUL terminators */
  return eslOK;
}

/* Function:  esl_msa_MinimGaps()
 * Synopsis:  Remove columns containing all gap symbols.
 * Incept:    SRE, Sun Feb 27 11:03:42 2005 [St. Louis]
 *
 * Purpose:   Remove all columns in the multiple alignment <msa>
 *            that consist entirely of gaps or missing data.
 *            
 *            For a text mode alignment, <gaps> is a string defining
 *            the gap characters, such as <"-_.">. For a digital mode
 *            alignment, <gaps> may be passed as <NULL>, because the
 *            internal alphabet already knows what the gap and missing
 *            data characters are.
 *            
 *            <msa> is changed in-place to a narrower alignment
 *            containing fewer columns. All per-residue and per-column
 *            annotation is altered appropriately for the columns that
 *            remain in the new alignment.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEMEM> on allocation failure.
 *            Possibilities from esl_msa_ColumnSubset() call:
 *            <eslESYNTAX> if WUSS string for SS_cons or msa->ss 
 *            following esl_wuss_nopseudo() is inconsistent.
 *            <eslEINVAL> if a derived ct array implies a pknotted SS.
 *
 * Xref:      squid's MSAMingap().
 */
int
esl_msa_MinimGaps(ESL_MSA *msa, char *errbuf, const char *gaps)
{
  int    *useme = NULL;	/* array of TRUE/FALSE flags for which cols to keep */
  int64_t apos;		/* column index   */
  int     idx;		/* sequence index */
  int     status;

  ESL_ALLOC(useme, sizeof(int) * (msa->alen+1)); /* +1 is just to deal w/ alen=0 special case */

#ifdef eslAUGMENT_ALPHABET	   /* digital mode case */
  if (msa->flags & eslMSA_DIGITAL) /* be careful of off-by-one: useme is 0..L-1 indexed */
    {
      for (apos = 1; apos <= msa->alen; apos++)
	{
	  for (idx = 0; idx < msa->nseq; idx++)
	    if (! esl_abc_XIsGap    (msa->abc, msa->ax[idx][apos]) &&
		! esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos]))
	      break;
	  if (idx == msa->nseq) useme[apos-1] = FALSE; else useme[apos-1] = TRUE;
	}
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL)) /* text mode case */
    {
      for (apos = 0; apos < msa->alen; apos++)
	{
	  for (idx = 0; idx < msa->nseq; idx++)
	    if (strchr(gaps, msa->aseq[idx][apos]) == NULL)
	      break;
	  if (idx == msa->nseq) useme[apos] = FALSE; else useme[apos] = TRUE;
	}
    }

  if((status = esl_msa_ColumnSubset(msa, errbuf, useme)) != eslOK) return status;
  free(useme);
  return eslOK;

 ERROR:
  if (useme != NULL) free(useme);
  return status;
}

/* Function:  esl_msa_NoGaps()
 * Synopsis:  Remove columns containing any gap symbol.
 * Incept:    SRE, Sun Feb 27 10:17:58 2005 [St. Louis]
 *
 * Purpose:   Remove all columns in the multiple alignment <msa> that
 *            contain any gaps or missing data, such that the modified
 *            MSA consists only of ungapped columns (a solid block of
 *            residues). 
 *            
 *            This is useful for filtering alignments prior to
 *            phylogenetic analysis using programs that can't deal
 *            with gaps.
 *            
 *            For a text mode alignment, <gaps> is a string defining
 *            the gap characters, such as <"-_.">. For a digital mode
 *            alignment, <gaps> may be passed as <NULL>, because the
 *            internal alphabet already knows what the gap and
 *            missing data characters are.
 *    
 *            <msa> is changed in-place to a narrower alignment
 *            containing fewer columns. All per-residue and per-column
 *            annotation is altered appropriately for the columns that
 *            remain in the new alignment.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            Possibilities from esl_msa_ColumnSubset() call:
 *            <eslESYNTAX> if WUSS string for SS_cons or msa->ss 
 *            following esl_wuss_nopseudo() is inconsistent.
 *            <eslEINVAL> if a derived ct array implies a pknotted SS.
 *
 * Xref:      squid's MSANogap().
 */
int
esl_msa_NoGaps(ESL_MSA *msa, char *errbuf, const char *gaps)
{
  int    *useme = NULL;	/* array of TRUE/FALSE flags for which cols to keep */
  int64_t apos;		/* column index */
  int     idx;		/* sequence index */
  int     status;

  ESL_ALLOC(useme, sizeof(int) * (msa->alen+1)); /* +1 is only to deal with alen=0 special case */

#ifdef eslAUGMENT_ALPHABET	   /* digital mode case */
  if (msa->flags & eslMSA_DIGITAL) /* be careful of off-by-one: useme is 0..L-1 indexed */
    {
      for (apos = 1; apos <= msa->alen; apos++)
	{
	  for (idx = 0; idx < msa->nseq; idx++)
	    if (esl_abc_XIsGap    (msa->abc, msa->ax[idx][apos]) ||
		esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos]))
	      break;
	  if (idx == msa->nseq) useme[apos-1] = TRUE; else useme[apos-1] = FALSE;
	}
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL)) /* text mode case */
    {
      for (apos = 0; apos < msa->alen; apos++)
	{
	  for (idx = 0; idx < msa->nseq; idx++)
	    if (strchr(gaps, msa->aseq[idx][apos]) != NULL)
	      break;
	  if (idx == msa->nseq) useme[apos] = TRUE; else useme[apos] = FALSE;
	}
    }

  esl_msa_ColumnSubset(msa, errbuf, useme);
  free(useme);
  return eslOK;

 ERROR:
  if (useme != NULL) free(useme);
  return status;
}


/* Function:  esl_msa_SymConvert()
 * Synopsis:  Global search/replace of symbols in an MSA.
 * Incept:    SRE, Sun Feb 27 11:20:41 2005 [St. Louis]
 *
 * Purpose:   In the aligned sequences in a text-mode <msa>, convert any
 *            residue in the string <oldsyms> to its counterpart (at the same
 *            position) in string <newsyms>.
 * 
 *            To convert DNA to RNA, <oldsyms> could be "Tt" and
 *            <newsyms> could be "Uu". To convert IUPAC symbols to
 *            N's, <oldsyms> could be "RYMKSWHBVDrymkswhbvd" and
 *            <newsyms> could be "NNNNNNNNNNnnnnnnnnnn". 
 *            
 *            As a special case, if <newsyms> consists of a single
 *            character, then any character in the <oldsyms> is 
 *            converted to this character. 
 *            
 *            Thus, <newsyms> must either be of the same length as
 *            <oldsyms>, or of length 1. Anything else will cause
 *            undefined behavior (and probably segfault). 
 *            
 *            The conversion is done in-place, so the <msa> is
 *            modified.
 *            
 *            This is a poor man's hack for processing text mode MSAs
 *            into a more consistent text alphabet. It is unnecessary
 *            for digital mode MSAs, which are already in a standard
 *            internal alphabet. Calling <esl_msa_SymConvert()> on a
 *            digital mode alignment throws an <eslEINVAL> error.
 *            
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEINVAL> if <msa> is in digital mode, or if the <oldsyms>
 *            and <newsyms> strings aren't valid together.
 */
int
esl_msa_SymConvert(ESL_MSA *msa, const char *oldsyms, const char *newsyms)
{
  int64_t apos;			/* column index */
  int     idx;			/* sequence index */
  char   *sptr;
  int     special;

  if (msa->flags & eslMSA_DIGITAL)
    ESL_EXCEPTION(eslEINVAL, "can't SymConvert on digital mode alignment");
  if ((strlen(oldsyms) != strlen(newsyms)) && strlen(newsyms) != 1)
    ESL_EXCEPTION(eslEINVAL, "invalid newsyms/oldsyms pair");

  special = (strlen(newsyms) == 1 ? TRUE : FALSE);

  for (apos = 0; apos < msa->alen; apos++)
    for (idx = 0; idx < msa->nseq; idx++)
      if ((sptr = strchr(oldsyms, msa->aseq[idx][apos])) != NULL)
	msa->aseq[idx][apos] = (special ? *newsyms : newsyms[sptr-oldsyms]);
  return eslOK;
}

/* Function:  esl_msa_AddComment()
 * Incept:    SRE, Tue Jun  1 17:37:21 1999 [St. Louis]
 *
 * Purpose:   Add an (unparsed) comment line to the MSA structure, 
 *            allocating as necessary.
 *
 * Args:      msa - a multiple alignment
 *            s   - comment line to add
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_msa_AddComment(ESL_MSA *msa, char *s)
{
  void *p;
  int   status;

  /* If this is our first recorded comment, we need to allocate;
   * and if we've filled available space, we need to reallocate.
   */
  if (msa->comment == NULL) {
    ESL_ALLOC(msa->comment, sizeof(char *) * 16);
    msa->alloc_ncomment = 16;
  }
  if (msa->ncomment == msa->alloc_ncomment) {
    ESL_RALLOC(msa->comment, p, sizeof(char *) * msa->alloc_ncomment * 2);
    msa->alloc_ncomment *= 2;
  }
  if ((status = esl_strdup(s, -1, &(msa->comment[msa->ncomment]))) != eslOK) goto ERROR;
  msa->ncomment++;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  esl_msa_AddGF()
 * Incept:    SRE, Tue Jun  1 17:37:21 1999 [St. Louis]
 *
 * Purpose:   Add an unparsed \verb+#=GF+ markup line to the MSA, 
 *            allocating as necessary. <tag> is the GF markup 
 *            tag; <value> is the text associated w/ that tag.
 *
 * Args:      msa - a multiple alignment
 *            tag - markup tag 
 *            value - markup text
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_msa_AddGF(ESL_MSA *msa, char *tag, char *value)
{  
  void *p;
  int   n;
  int   status;

  /* If this is our first recorded unparsed #=GF line, we need to allocate().
   */
  if (msa->gf_tag == NULL) {
    ESL_ALLOC(msa->gf_tag, sizeof(char *) * 16);
    ESL_ALLOC(msa->gf,     sizeof(char *) * 16);
    msa->alloc_ngf = 16;
  }
  /* or if we're out of room for new GF's, reallocate by doubling
   */
  if (msa->ngf == msa->alloc_ngf) {
    n = msa->alloc_ngf * 2;
    ESL_RALLOC(msa->gf_tag, p, sizeof(char *) * n);
    ESL_RALLOC(msa->gf,     p, sizeof(char *) * n);
    msa->alloc_ngf = n;
  }

  if ((status = esl_strdup(tag,  -1,  &(msa->gf_tag[msa->ngf]))) != eslOK) goto ERROR;
  if ((status = esl_strdup(value, -1, &(msa->gf[msa->ngf])))     != eslOK) goto ERROR;
  msa->ngf++;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  esl_msa_AddGS()
 * Incept:    SRE, Tue Jun  1 17:37:21 1999 [St. Louis]
 *
 * Purpose:   Add an unparsed \verb+#=GS+ markup line to the MSA, 
 *            allocating as necessary. It's possible that we 
 *            could get more than one of the same type of GS 
 *            tag per sequence; for example, "DR PDB;" structure 
 *            links in Pfam.  Hack: handle these by appending to 
 *            the string, in a \verb+\n+ separated fashion.
 *
 * Args:      msa    - multiple alignment structure
 *            tag    - markup tag (e.g. "AC")
 *            sqidx  - index of sequence to assoc markup with (0..nseq-1)
 *            value  - markup (e.g. "P00666")
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_msa_AddGS(ESL_MSA *msa, char *tag, int sqidx, char *value)
{
  void *p;
  int   tagidx;
  int   i;
  int   status;

  /* first GS tag? init&allocate  */
  if (msa->gs_tag == NULL)	
    {
#ifdef eslAUGMENT_KEYHASH
      msa->gs_idx = esl_keyhash_Create();
      status = esl_key_Store(msa->gs_idx, tag, &tagidx);
      if (status != eslOK && status != eslEDUP) return status;
      ESL_DASSERT1((tagidx == 0));
#else
      tagidx = 0;
#endif
      ESL_ALLOC(msa->gs_tag, sizeof(char *));  /* one at a time. */
      ESL_ALLOC(msa->gs,     sizeof(char **));
      ESL_ALLOC(msa->gs[0],  sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++)
	msa->gs[0][i] = NULL;
    }
  else 
    {
      /* Get a tagidx for this GS tag.
       * tagidx < ngs; we already saw this tag;
       * tagidx == ngs; this is a new one.
       */
#ifdef eslAUGMENT_KEYHASH
      status = esl_key_Store(msa->gs_idx, tag, &tagidx);
      if (status != eslOK && status != eslEDUP) return status;
#else
      for (tagidx = 0; tagidx < msa->ngs; tagidx++)
	if (strcmp(msa->gs_tag[tagidx], tag) == 0) break;
#endif
      /* Reallocation (in blocks of 1) */
      if (tagidx == msa->ngs ) 
	{
	  ESL_RALLOC(msa->gs_tag, p, (msa->ngs+1) * sizeof(char *));
	  ESL_RALLOC(msa->gs,     p, (msa->ngs+1) * sizeof(char **));
	  ESL_ALLOC(msa->gs[msa->ngs], sizeof(char *) * msa->sqalloc);
	  for (i = 0; i < msa->sqalloc; i++) 
	    msa->gs[msa->ngs][i] = NULL;
	}
    }

  /* Store the tag, if it's new.
   */
  if (tagidx == msa->ngs) 
    {
      if ((status = esl_strdup(tag, -1, &(msa->gs_tag[tagidx]))) != eslOK) goto ERROR;
      msa->ngs++;
    }
  
  /* Store the annotation on the sequence.
   * If seq is unannotated, dup the value; if
   * seq already has a GS annotation, cat a \n, then cat the value.
   */
  if (msa->gs[tagidx][sqidx] == NULL)
    {
      if ((status = esl_strdup(value, -1, &(msa->gs[tagidx][sqidx]))) != eslOK) goto ERROR;
    }
  else 
    {			
      int n1,n2;
      n1 = strlen(msa->gs[tagidx][sqidx]);
      n2 = strlen(value);
      ESL_RALLOC(msa->gs[tagidx][sqidx], p, sizeof(char) * (n1+n2+2));
      msa->gs[tagidx][sqidx][n1] = '\n';
      strcpy(msa->gs[tagidx][sqidx]+n1+1, value);
    }
  return eslOK;

 ERROR:
  return status;
} 

/* Function:  esl_msa_AppendGC()
 * Incept:    SRE, Tue Jun  1 17:37:21 1999 [St. Louis]
 *
 * Purpose:   Add an unparsed \verb+#=GC+ markup line to the MSA 
 *            structure, allocating as necessary. When called 
 *            multiple times for the same tag, appends value 
 *            strings together -- used when parsing multiblock 
 *            alignment files, for example.
 *
 * Args:      msa   - multiple alignment structure
 *            tag   - markup tag (e.g. "CS")
 *            value - markup, one char per aligned column      
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_msa_AppendGC(ESL_MSA *msa, char *tag, char *value)
{
  int   tagidx;
  int   status;
  void *p;

  /* Is this an unparsed tag name that we recognize?
   * If not, handle adding it to index, and reallocating
   * as needed.
   */
  if (msa->gc_tag == NULL)	/* first tag? init&allocate  */
    {
#ifdef eslAUGMENT_KEYHASH
      msa->gc_idx = esl_keyhash_Create();
      status = esl_key_Store(msa->gc_idx, tag, &tagidx);      
      if (status != eslOK && status != eslEDUP) return status;
      ESL_DASSERT1((tagidx == 0));
#else
      tagidx = 0;
#endif
      ESL_ALLOC(msa->gc_tag, sizeof(char **));
      ESL_ALLOC(msa->gc,     sizeof(char **));
      msa->gc[0]  = NULL;
    }
  else
    {			/* new tag? */
      /* get tagidx for this GC tag. existing tag: <ngc; new: == ngc. */
#ifdef eslAUGMENT_KEYHASH
      status = esl_key_Store(msa->gc_idx, tag, &tagidx);
      if (status != eslOK && status != eslEDUP) goto ERROR;
#else
      for (tagidx = 0; tagidx < msa->ngc; tagidx++)
	if (strcmp(msa->gc_tag[tagidx], tag) == 0) break;
#endif
      /* Reallocate, in block of one tag at a time
       */
      if (tagidx == msa->ngc)
	{
	  ESL_RALLOC(msa->gc_tag, p, (msa->ngc+1) * sizeof(char **));
	  ESL_RALLOC(msa->gc,     p, (msa->ngc+1) * sizeof(char **));
	  msa->gc[tagidx] = NULL;
	}
    }
  /* new tag? store it.
   */
  if (tagidx == msa->ngc) 
    {
      if ((status = esl_strdup(tag, -1, &(msa->gc_tag[tagidx]))) != eslOK) goto ERROR;
      msa->ngc++;
    }
  return (esl_strcat(&(msa->gc[tagidx]), -1, value, -1));

 ERROR:
  return status;
}

/* Function:  esl_msa_AppendGR()
 * Incept:    SRE, Thu Jun  3 06:34:38 1999 [Madison]
 *
 * Purpose:   Add an unparsed \verb+#=GR+ markup line to the MSA structure, 
 *            allocating as necessary.
 *              
 *            When called multiple times for the same tag, appends 
 *            value strings together -- used when parsing multiblock 
 *            alignment files, for example.
 *
 * Args:      msa    - multiple alignment structure
 *            tag    - markup tag (e.g. "SS")
 *            sqidx  - index of seq to assoc markup with (0..nseq-1)
 *            value  - markup, one char per aligned column      
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_msa_AppendGR(ESL_MSA *msa, char *tag, int sqidx, char *value)
{
  void *p;
  int tagidx;
  int i;
  int status;

  if (msa->gr_tag == NULL)	/* first tag? init&allocate  */
    {
#ifdef eslAUGMENT_KEYHASH
      msa->gr_idx = esl_keyhash_Create();
      status = esl_key_Store(msa->gr_idx, tag, &tagidx);
      if (status != eslOK && status != eslEDUP) return status;
      ESL_DASSERT1((tagidx == 0));
#else
      tagidx = 0;
#endif
      ESL_ALLOC(msa->gr_tag, sizeof(char *));
      ESL_ALLOC(msa->gr,     sizeof(char **));
      ESL_ALLOC(msa->gr[0],  sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++) 
	msa->gr[0][i] = NULL;
    }
  else 
    {
      /* get tagidx for this GR tag. existing<ngr; new=ngr.
       */
#ifdef eslAUGMENT_KEYHASH
      status = esl_key_Store(msa->gr_idx, tag, &tagidx);
      if (status != eslOK && status != eslEDUP) return status;
#else
      for (tagidx = 0; tagidx < msa->ngr; tagidx++)
	if (strcmp(msa->gr_tag[tagidx], tag) == 0) break;
#endif
      /* if a new tag, realloc for it */      
      if (tagidx == msa->ngr)
	{ 
	  ESL_RALLOC(msa->gr_tag, p, (msa->ngr+1) * sizeof(char *));
	  ESL_RALLOC(msa->gr,     p, (msa->ngr+1) * sizeof(char **));
	  ESL_ALLOC(msa->gr[msa->ngr], sizeof(char *) * msa->sqalloc);
	  for (i = 0; i < msa->sqalloc; i++) 
	    msa->gr[msa->ngr][i] = NULL;
	}
    }

  if (tagidx == msa->ngr) 
    {
      if ((status = esl_strdup(tag, -1, &(msa->gr_tag[tagidx]))) != eslOK) goto ERROR;
      msa->ngr++;
    }
  return (esl_strcat(&(msa->gr[tagidx][sqidx]), -1, value, -1));

 ERROR:
  return status;
}


/* Function:  esl_msa_Checksum()
 * Synopsis:  Calculate a checksum for an MSA.
 * Incept:    SRE, Tue Sep 16 13:23:34 2008 [Janelia]
 *
 * Purpose:   Calculates a 32-bit checksum for <msa>.
 * 
 *            Only the alignment data are considered, not the sequence
 *            names or other annotation. For text mode alignments, the
 *            checksum is case sensitive.
 *            
 *            This is used as a quick way to try to verify that a
 *            given alignment is identical to an expected one; for
 *            example, when HMMER is mapping new sequence alignments
 *            onto exactly the same seed alignment an HMM was built
 *            from.
 *
 * Returns:   <eslOK> on success.
 *
 * Xref:      The checksum is a modified version of Jenkin's hash;
 *            see <esl_keyhash> for the original and citations.
 */
int
esl_msa_Checksum(const ESL_MSA *msa, uint32_t *ret_checksum)
{
  uint32_t val = 0;
  int      i,pos;

#ifdef eslAUGMENT_ALPHABET
  if (msa->flags & eslMSA_DIGITAL)
    {
      for (i = 0; i < msa->nseq; i++)
	for (pos = 0; pos <= msa->alen; pos++)
	  {
	    val += msa->ax[i][pos];
	    val += (val << 10);
	    val ^= (val >>  6);
	  }
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL))
    {
      for (i = 0; i < msa->nseq; i++)
	for (pos = 0; pos <= msa->alen; pos++)
	  {
	    val += msa->ax[i][pos];
	    val += (val << 10);
	    val ^= (val >>  6);
	  }
    }
  val += (val <<  3);
  val ^= (val >> 11);
  val += (val << 15);

  *ret_checksum = val;
  return eslOK;
}
/*-------------------- end of misc MSA functions ----------------------*/





/******************************************************************************
 * 7. Debugging/development routines.
 *****************************************************************************/

/* Function:  esl_msa_CreateFromString()
 * Synopsis:  Creates a small <ESL_MSA> from a test case string.
 * Incept:    SRE, Sat Nov 11 12:09:04 2006 [Janelia]
 *
 * Purpose:   A convenience for making small test cases in the test
 *            suites: given the contents of a complete multiple
 *            sequence alignment file as a single string <s> in
 *            alignment format <fmt>, convert it to an <ESL_MSA>.
 *            
 *            For example, 
 *            {\small\begin{verbatim}
 *            esl_msa_CreateFromString("# STOCKHOLM 1.0\n\nseq1 AAAAA\nseq2 AAAAA\n//\n", 
 *                                     eslMSAFILE_STOCKHOLM)
 *            \end{verbatim}}
 *            creates an ungapped alignment of two AAAAA sequences.
 *
 * Returns:   a pointer to the new <ESL_MSA> on success.
 *
 * Throws:    <NULL> if it fails to obtain, open, or read the temporary file
 *            that it puts the string <s> in.
 */
ESL_MSA *
esl_msa_CreateFromString(const char *s, int fmt)
{
  char         tmpfile[16] = "esltmpXXXXXX";
  FILE        *fp          = NULL;
  ESL_MSAFILE *mfp         = NULL;
  ESL_MSA     *msa         = NULL;

  if (esl_tmpfile_named(tmpfile, &fp)            != eslOK) goto ERROR;
  fprintf(fp, s);
  fclose(fp); 
  fp = NULL;
  if (esl_msafile_Open(tmpfile, fmt, NULL, &mfp) != eslOK) goto ERROR;
  if (esl_msa_Read(mfp, &msa)                    != eslOK) goto ERROR;

  esl_msafile_Close(mfp);
  remove(tmpfile);
  return msa;

 ERROR:
  if (fp  != NULL) fclose(fp);
  if (mfp != NULL) esl_msafile_Close(mfp);
  if (strcmp(tmpfile, "esltmpXXXXXX") != 0) remove(tmpfile);
  if (msa != NULL) esl_msa_Destroy(msa);                        
  return NULL;
}

/* Function:  esl_msa_Compare()
 * Synopsis:  Compare two MSAs for equality.
 * Incept:    SRE, Wed Jun 13 10:40:05 2007 [Janelia]
 *
 * Purpose:   Returns <eslOK> if the mandatory and optional contents
 *            of MSAs <a1> and <a2> are identical; otherwise return
 *            <eslFAIL>.
 *            
 *            Only mandatory and parsed optional information is
 *            compared. Unparsed Stockholm markup is not compared.
 */
int
esl_msa_Compare(ESL_MSA *a1, ESL_MSA *a2)
{
  if (esl_msa_CompareMandatory(a1, a2) != eslOK) return eslFAIL;
  if (esl_msa_CompareOptional(a1, a2)  != eslOK) return eslFAIL;
  return eslOK;
}

/* Function:  esl_msa_CompareMandatory()
 * Synopsis:  Compare mandatory subset of MSA contents.
 * Incept:    SRE, Wed Jun 13 09:42:56 2007 [Janelia]
 *
 * Purpose:   Compare mandatory contents of two MSAs, <a1> and <a2>.
 *            This comprises <aseq> (or <ax>, for a digital alignment);
 *            <sqname>, <wgt>, <alen>, <nseq>, and <flags>.
 *
 * Returns:   <eslOK> if the MSAs are identical; 
 *            <eslFAIL> if they are not.
 */
int
esl_msa_CompareMandatory(ESL_MSA *a1, ESL_MSA *a2)
{
  int i;

  if (a1->nseq  != a2->nseq)  return eslFAIL;
  if (a1->alen  != a2->alen)  return eslFAIL;
  if (a1->flags != a2->flags) return eslFAIL;

  for (i = 0; i < a1->nseq; i++)
    {
      if (strcmp(a1->sqname[i], a2->sqname[i])        != 0)     return eslFAIL;
      if (esl_DCompare(a1->wgt[i], a2->wgt[i], 0.001) != eslOK) return eslFAIL;
#ifdef eslAUGMENT_ALPHABET
      if ((a1->flags & eslMSA_DIGITAL) &&
	  memcmp(a1->ax[i], a2->ax[i], sizeof(ESL_DSQ) * (a1->alen+2)) != 0) 
	return eslFAIL;
#endif
      if (! (a1->flags & eslMSA_DIGITAL) && strcmp(a1->aseq[i], a2->aseq[i]) != 0) return eslFAIL;
    }
  return eslOK;
}

/* Function:  esl_msa_CompareOptional()
 * Synopsis:  Compare optional subset of MSA contents.
 * Incept:    SRE, Wed Jun 13 09:52:48 2007 [Janelia]
 *
 * Purpose:   Compare optional contents of two MSAs, <a1> and <a2>.
 *
 * Returns:   <eslOK> if the MSAs are identical; 
 *            <eslFAIL> if they are not.
 */
int
esl_msa_CompareOptional(ESL_MSA *a1, ESL_MSA *a2)
{
  int i;

  if (esl_CCompare(a1->name,    a2->name)    != eslOK) return eslFAIL;
  if (esl_CCompare(a1->desc,    a2->desc)    != eslOK) return eslFAIL;
  if (esl_CCompare(a1->acc,     a2->acc)     != eslOK) return eslFAIL;
  if (esl_CCompare(a1->au,      a2->au)      != eslOK) return eslFAIL;
  if (esl_CCompare(a1->ss_cons, a2->ss_cons) != eslOK) return eslFAIL;
  if (esl_CCompare(a1->sa_cons, a2->sa_cons) != eslOK) return eslFAIL;
  if (esl_CCompare(a1->rf,      a2->rf)      != eslOK) return eslFAIL;
  
  if (a1->sqacc != NULL && a2->sqacc != NULL) {
    for (i = 0; i < a1->nseq; i++) if (esl_CCompare(a1->sqacc[i], a2->sqacc[i]) != eslOK) return eslFAIL;
  } else if (a1->sqacc != NULL || a2->sqacc != NULL) return eslFAIL;

  if (a1->sqdesc != NULL && a2->sqdesc != NULL) {
    for (i = 0; i < a1->nseq; i++) if (esl_CCompare(a1->sqdesc[i], a2->sqdesc[i]) != eslOK) return eslFAIL;
  } else if (a1->sqdesc != NULL || a2->sqdesc != NULL) return eslFAIL;

  if (a1->ss != NULL && a2->ss != NULL) {
    for (i = 0; i < a1->nseq; i++) if (esl_CCompare(a1->ss[i], a2->ss[i]) != eslOK) return eslFAIL;
  } else if (a1->ss != NULL || a2->ss != NULL) return eslFAIL;

  if (a1->sa != NULL && a2->sa != NULL) {
    for (i = 0; i < a1->nseq; i++) if (esl_CCompare(a1->sa[i], a2->sa[i]) != eslOK) return eslFAIL;
  } else if (a1->sa != NULL || a2->sa != NULL) return eslFAIL;
  
  for (i = 0; i < eslMSA_NCUTS; i++)
    {
      if (a1->cutset[i] && a2->cutset[i]) {
	if (esl_FCompare(a1->cutoff[i], a2->cutoff[i], 0.01) != eslOK) return eslFAIL;
      } else if (a1->cutset[i] || a2->cutset[i]) return eslFAIL;
    }
  return eslOK;
}


/*---------------- end of debugging/development routines  -------------------*/



/******************************************************************************
 * 8. Benchmark driver.
 *****************************************************************************/
#ifdef eslMSA_BENCHMARK
/* gcc -O2 -o msa_benchmark -I. -L. -DeslMSA_BENCHMARK esl_msa.c -leasel -lm
 * ./benchmark Pfam
 */
#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_stopwatch.h"

static ESL_OPTIONS options[] = {
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-v",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "be verbose: show info as each msa is read",        0 },
  { 0,0,0,0,0,0,0,0 },
};
static char usage[]  = "[-options] <msafile>";
static char banner[] = "benchmarking speed of MSA reading";

int
main(int argc, char **argv)
{
  ESL_GETOPTS   *go        = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_STOPWATCH *w         = esl_stopwatch_Create();
  char          *filename  = esl_opt_GetArg(go, 1);
  int            fmt       = eslMSAFILE_UNKNOWN;
  ESL_MSAFILE   *afp       = NULL;
  ESL_MSA       *msa       = NULL;
  int            status;
  int            nali;
  int            alphatype;

  /* Open the alignment file in text mode */
  status = esl_msafile_Open(filename, fmt, NULL, &afp);
  if (status == eslENOTFOUND)     esl_fatal("Alignment file %s doesn't exist or is not readable\n", filename);
  else if (status == eslEFORMAT)  esl_fatal("Couldn't determine format of alignment %s\n", filename);
  else if (status != eslOK)       esl_fatal("Alignment file open failed with error %d\n", status);

  /* Loop over all the alignments. */
  esl_stopwatch_Start(w);
  nali = 0;
  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;

      if (esl_opt_GetBoolean(go, "-v")) 
	{
	  esl_msa_GuessAlphabet(msa, &alphatype);
	  printf("%5d %15s %6d %5d %7s\n", 
		 nali, msa->name, msa->nseq, msa->alen, esl_abc_DescribeType(alphatype));
	}

      esl_msa_Destroy(msa);
    }
  if (status == eslEFORMAT)
    esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", 
	      afp->linenumber, afp->fname, afp->errbuf, afp->buf);
  else if (status != eslEOF)
    esl_fatal("Alignment file read failed with error code %d\n", status);

  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "CPU Time: ");

  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  esl_msafile_Close(afp);
  exit(0);
}
#endif /*eslMSA_BENCHMARK*/
/*---------------------- end of benchmark driver ----------------------*/



/******************************************************************************
 * 9. Unit tests
 *****************************************************************************/
#ifdef eslMSA_TESTDRIVE

/* write_known_msa()
 * Write a known MSA to a tmpfile in Stockholm format.
 */
static void
write_known_msa(FILE *ofp)
{
  fprintf(ofp, "# STOCKHOLM 1.0\n");
  fprintf(ofp, "seq1 --ACDEFGHIK~LMNPQRS-TVWY\n");
  fprintf(ofp, "seq2 aaACDEFGHIK~LMNPQRS-TVWY\n");
  fprintf(ofp, "seq3 aaACDEFGHIK~LMNPQRS-TVWY\n");
  fprintf(ofp, "\n");
  fprintf(ofp, "seq1 ACDEFGHIKLMNPQRSTVWY~~~\n");
  fprintf(ofp, "seq2 ACDEFGHIKLMNPQRSTVWYyyy\n");
  fprintf(ofp, "seq3 ACDEFGHIKLMNPQRSTVWYyyy\n");
  fprintf(ofp, "//\n");
  return;
}
  
/* compare_to_known() 
 * SRE, Thu Sep  7 09:52:07 2006 [Janelia]
 * Spotcheck an ESL_MSA to make sure it matches the test known alignment.
 */
static void
compare_to_known(ESL_MSA *msa)
{
  if (msa->alen != 47)                     esl_fatal("bad alen");
  if (msa->nseq != 3)                      esl_fatal("bad nseq");
  if (strcmp(msa->sqname[1], "seq2") != 0) esl_fatal("bad sqname");
#ifdef eslAUGMENT_ALPHABET
  if (msa->flags & eslMSA_DIGITAL)
    {
      if (! esl_abc_XIsGap(msa->abc, msa->ax[0][2]))      esl_fatal("no gap where expected");
      if (! esl_abc_XIsMissing(msa->abc, msa->ax[0][47])) esl_fatal("no missing-data symbol where expected");
      if (msa->ax[1][1]  != 0)                            esl_fatal("spotcheck on ax failed"); /* 0=A */
      if (msa->ax[1][47] != 19)                           esl_fatal("spotcheck on ax failed"); /*19=Y */
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL))
    {
      if (strcasecmp(msa->aseq[0], "--ACDEFGHIK~LMNPQRS-TVWYACDEFGHIKLMNPQRSTVWY~~~") != 0) esl_fatal("aseq 0 is bad");
      if (strcasecmp(msa->aseq[1], "aaACDEFGHIK~LMNPQRS-TVWYACDEFGHIKLMNPQRSTVWYyyy") != 0) esl_fatal("aseq 1 is bad");
      if (strcasecmp(msa->aseq[2], "aaACDEFGHIK~LMNPQRS-TVWYACDEFGHIKLMNPQRSTVWYyyy") != 0) esl_fatal("aseq 2 is bad");
    }
  return;
}

/* msa_compare()
 * SRE, Fri Sep  8 08:13:38 2006 [Janelia]
 * 
 * Compares two MSAs; returns eslOK if they appear to be the same, eslFAIL if not.
 * (Not worth putting in external API. Only useful for testing purposes.)
 * Not a complete comparison: just checks mode, sequence names, and aligned data.
 * MSAs have to be in same mode (text vs. digital).
 */
static void
msa_compare(ESL_MSA *m1, ESL_MSA *m2)
{
  int i;

  if (m1->nseq  != m2->nseq   ||
      m1->alen  != m2->alen   ||
      m1->flags != m2->flags)
    esl_fatal("msa1, msa2 differ in nseq, alen, or flags");

  for (i = 0; i < m1->nseq; i++)
    {
      if (strcmp(m1->sqname[i], m2->sqname[i]) != 0) esl_fatal("msa1, msa2 sqnames differ for %d", i);
#ifdef eslAUGMENT_ALPHABET
      if ((m1->flags & eslMSA_DIGITAL) && 
	  memcmp(m1->ax[i], m2->ax[i], sizeof(ESL_DSQ) * (m1->alen+2)) != 0) 
	esl_fatal("msa1, msa2 digital sequences differ for %d", i);
#endif
      if (! (m1->flags & eslMSA_DIGITAL) && 
	  strcmp(m1->aseq[i], m2->aseq[i]) != 0) 
	esl_fatal("msa1, msa2 sequences differ for %d", i);
    }
  return;
}

/* Unit tests for every function in the exposed API
 */
static void
utest_Create(void)
{
  ESL_MSA *msa = NULL;

  msa = esl_msa_Create(16, -1);	  /* nseq blocksize 16, growable */
  esl_msa_Destroy(msa);
  msa = esl_msa_Create(16, 100);  /* nseq=16, alen=100, not growable */
  esl_msa_Destroy(msa);

  return;
}

static void
utest_Destroy(void)
{
  ESL_MSA *msa = NULL;
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET *abc;
#endif

  msa = esl_msa_Create(16, -1);	
  esl_msa_Destroy(msa);	 	  /* normal usage */

#ifdef eslAUGMENT_ALPHABET
  abc = esl_alphabet_Create(eslRNA);
  msa = esl_msa_CreateDigital(abc, 16, 100);	
  esl_msa_Destroy(msa);	 	  /* normal usage, digital mode */
  esl_alphabet_Destroy(abc);
#endif

  esl_msa_Destroy(NULL);	  /* should tolerate NULL argument */
  return;
}

static void
utest_Expand(void)
{
  ESL_MSA *msa = NULL;
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET *abc;
#endif

  msa = esl_msa_Create(16, -1);                	    /* growable */
  if (esl_msa_Expand(msa) != eslOK) esl_fatal("Expand failed"); /* expand by 2x in nseq */
  esl_msa_Destroy(msa);

  msa = esl_msa_Create(16, 100);                        /* not growable */
#ifdef eslTEST_THROWING
  if (esl_msa_Expand(msa) != eslEINVAL) esl_fatal("Expand should have failed but didn't"); /* should fail w/ EINVAL*/
#endif
  esl_msa_Destroy(msa);
  
#ifdef eslAUGMENT_ALPHABET
  abc = esl_alphabet_Create(eslDNA);
  msa = esl_msa_CreateDigital(abc, 16, -1);               /* growable */
  if (esl_msa_Expand(msa) != eslOK) esl_fatal("Expand failed"); /* expand by 2x in nseq */
  esl_msa_Destroy(msa);

  msa = esl_msa_CreateDigital(abc, 16, 100);                 /* not growable */
#ifdef eslTEST_THROWING
  if (esl_msa_Expand(msa) != eslEINVAL) esl_fatal("Expand should have failed but didn't"); /* should fail w/ EINVAL*/
#endif /* eslTEST_THROWING*/
  esl_msa_Destroy(msa);
  esl_alphabet_Destroy(abc);
#endif
  return;
}

static void
utest_Open(char *tmpfile)	/* filename must be in /tmp */
{
  char        *msg      = "Open() unit test failed";
  ESL_MSAFILE *msafp    = NULL;
  int          status;
  
  status = esl_msafile_Open(tmpfile, eslMSAFILE_UNKNOWN, NULL, &msafp); 
  if (status != eslOK) esl_fatal(msg);
  esl_msafile_Close(msafp);

  status = esl_msafile_Open(tmpfile, eslMSAFILE_STOCKHOLM, NULL, &msafp);
  if (status != eslOK) esl_fatal(msg);
  esl_msafile_Close(msafp);

#ifdef HAVE_PUTENV
  putenv("ESLTEST=./");
  esl_FileTail(tmpfile, FALSE, &filename);
  status = esl_msafile_Open(filename, eslMSAFILE_STOCKHOLM, "ESLTEST", &msafp);
  if (status != eslOK) esl_fatal(msg);
  esl_msafile_Close(msafp);
  free(filename);
#endif

  return;
}

static void
utest_Close(char *filename)
{
  ESL_MSAFILE *msafp    = NULL;
  int status;

  status = esl_msafile_Open(filename, eslMSAFILE_UNKNOWN, NULL, &msafp); 
  if (status != eslOK) esl_fatal("Close() unit test failed");
  esl_msafile_Close(msafp);
  esl_msafile_Close(NULL);	/* should work */
  return;
}

#ifdef eslAUGMENT_ALPHABET
static void
utest_CreateDigital(ESL_ALPHABET *abc)
{
  char    *msg = "CreateDigital() unit test failure";
  ESL_MSA *msa = NULL;

  msa = esl_msa_CreateDigital(abc, 16, -1);	  /* nseq blocksize 16, growable */
  if (! (msa->flags & eslMSA_DIGITAL)) esl_fatal(msg);
  if (msa->ax   == NULL)               esl_fatal(msg);
  if (msa->aseq != NULL)               esl_fatal(msg);
  if (esl_msa_Expand(msa) != eslOK)    esl_fatal(msg);
  esl_msa_Destroy(msa);

  msa = esl_msa_CreateDigital(abc, 16, 100);  /* nseq=16, alen=100, not growable */
#ifdef eslTEST_THROWING
  if (esl_msa_Expand(msa) != eslEINVAL) esl_fatal(msg); /* shouldn't grow */
#endif
  esl_msa_Destroy(msa);

  return;
}
#endif /*eslAUGMENT_ALPHABET*/

#ifdef eslAUGMENT_ALPHABET
static void
utest_Digitize(ESL_ALPHABET *abc, char *filename)
{
  char        *msg = "Digitize() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  ESL_MSA     *msa = NULL;
  int c, i, pos;

  /* Get ourselves a copy of the known alignment that we can muck with */
  if (esl_msafile_Open(filename, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK)  esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK)                                       esl_fatal(msg);
  esl_msafile_Close(mfp);
  
  /* Deliberately corrupt it with inval character in the middle */
  i   = msa->nseq / 2;
  pos = msa->alen / 2;
  c   = msa->aseq[i][pos];
  msa->aseq[i][pos] = '%';
  if (esl_msa_Digitize(abc, msa) != eslEINVAL) esl_fatal(msg); /* should detect corruption as normal error */
  msa->aseq[i][pos] = c;	                               /* restore original         */
  compare_to_known(msa);
  if (esl_msa_Digitize(abc, msa) != eslOK)     esl_fatal(msg); /* should be fine now       */
  compare_to_known(msa);

  esl_msa_Destroy(msa);
  return;
}
#endif /*eslAUGMENT_ALPHABET*/


#ifdef eslAUGMENT_ALPHABET
static void
utest_Textize(ESL_ALPHABET *abc, char *filename)
{
  char        *msg = "Textize() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  ESL_MSA     *msa = NULL;

  if (esl_msafile_OpenDigital(abc, filename, eslMSAFILE_UNKNOWN, NULL, &mfp) != eslOK)  esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK)   esl_fatal(msg);
  if (esl_msa_Textize(msa)    != eslOK)   esl_fatal(msg);
  compare_to_known(msa);

  esl_msafile_Close(mfp);
  esl_msa_Destroy(msa);
  return;
}
#endif /*eslAUGMENT_ALPHABET*/

#ifdef eslAUGMENT_ALPHABET
static void
utest_OpenDigital(ESL_ALPHABET *abc, char *filename)  /* filename must be in /tmp */
{
  char        *msg   = "OpenDigital() unit test failure";
  ESL_MSAFILE *msafp = NULL;
  
  if (esl_msafile_OpenDigital(abc, filename, eslMSAFILE_UNKNOWN,   NULL, &msafp) != eslOK) esl_fatal(msg);  esl_msafile_Close(msafp);
  if (esl_msafile_OpenDigital(abc, filename, eslMSAFILE_STOCKHOLM, NULL, &msafp) != eslOK) esl_fatal(msg);  esl_msafile_Close(msafp);
#ifdef HAVE_PUTENV
  putenv("ESLTEST=./");
  esl_FileTail(tmpfile, FALSE, &filename);
  if (esl_msafile_OpenDigital(abc, filename, eslMSAFILE_STOCKHOLM, "ESLTEST", &msafp) != eslOK) esl_fatal(msg);
  esl_msafile_Close(msafp);
  free(filename);
#endif
  return;
}
#endif /*eslAUGMENT_ALPHABET*/

static void
utest_Read(char *filename)
{
  char        *msg = "Read() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  ESL_MSA     *msa = NULL;

  if (esl_msafile_Open(filename, eslMSAFILE_UNKNOWN, NULL, &mfp) != eslOK)  esl_fatal(msg);  
  if (esl_msa_Read(mfp, &msa) != eslOK)  esl_fatal(msg);
  compare_to_known(msa);
  esl_msa_Destroy(msa);

  if (esl_msa_Read(mfp, &msa) != eslEOF) esl_fatal(msg);
  if (msa != NULL)                       esl_fatal(msg);

  esl_msafile_Close(mfp);
  return;
}

static void
utest_Write(ESL_MSA *msa1)
{
  char        *msg  = "Write() unit test failure";
  ESL_MSAFILE *mfp  = NULL;
  ESL_MSA     *msa2 = NULL;
  FILE        *fp   = NULL;
  int      i;
  int      formats[] = { eslMSAFILE_STOCKHOLM, eslMSAFILE_PFAM, -1 }; /* -1=sentinel */
  char     template[16] = "esltmpXXXXXX";
  char     tmpfile[16];

  for (i = 0; formats[i] != -1; i++)
    {
      strcpy(tmpfile, template);
      if (esl_tmpfile_named(tmpfile, &fp) != eslOK)     esl_fatal(msg);
      if (esl_msa_Write(fp, msa1, formats[i]) != eslOK) esl_fatal(msg);
      fclose(fp);
  
#ifdef eslAUGMENT_ALPHABET
      if ((msa1->flags & eslMSA_DIGITAL) &&
	  esl_msafile_OpenDigital(msa1->abc, tmpfile, formats[i], NULL, &mfp) != eslOK)	esl_fatal(msg);
#endif
      if (! (msa1->flags & eslMSA_DIGITAL) &&
	  esl_msafile_Open(tmpfile, formats[i], NULL, &mfp) != eslOK) esl_fatal(msg);

      if (esl_msa_Read(mfp, &msa2) != eslOK) esl_fatal(msg);
      msa_compare(msa1, msa2);
      
      esl_msafile_Close(mfp);
      esl_msa_Destroy(msa2);
      remove(tmpfile);      
    }      
  return;
}

static void
utest_GuessFileFormat(void)
{
  /* SRE: To be filled in. Currently, esl_msa_GuessFileFormat() is a placeholder that
   * always guesses Stockholm
   */
  return;
}


static void
utest_SequenceSubset(ESL_MSA *m1)
{
  char    *msg   = "SequenceSubset() unit test failure";
  ESL_MSA *m2    = NULL;
  int     *useme = NULL;
  int      i,j;
  int      n2;

  /* Make every other sequence (1,3..) get excluded from the subset */
  useme = malloc(m1->nseq * sizeof(int));
  for (i = 0, n2 = 0; i < m1->nseq; i++)
    if (i%2 == 0) { useme[i] = TRUE; n2++; }
    else          useme[i] = FALSE;

  if (esl_msa_SequenceSubset(m1, useme, &m2) != eslOK) esl_fatal(msg);
  if (m2->nseq != n2) esl_fatal(msg);
  
  for (i = 0, j = 0; i < m1->nseq; i++)
    {
      if (useme[i])
	{
	  if (strcmp(m1->sqname[i], m2->sqname[j]) != 0) esl_fatal(msg);
	  if (! (m1->flags & eslMSA_DIGITAL) && (strcmp(m1->aseq[i],   m2->aseq[j])  != 0)) esl_fatal(msg);
#ifdef eslAUGMENT_ALPHABET
	  if (  (m1->flags & eslMSA_DIGITAL) && memcmp(m1->ax[i], m2->ax[j], sizeof(ESL_DSQ) * (m1->alen+2)) != 0) esl_fatal(msg);
#endif
	  j++;
	}
    }  
  esl_msa_Destroy(m2);
  free(useme);
  return;
}

static void
utest_MinimGaps(char *tmpfile)
{
  char        *msg = "MinimGaps() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  ESL_MSA     *msa = NULL;
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET *abc = NULL;
#endif

  if (esl_msafile_Open(tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK)                                     esl_fatal(msg);
  esl_msafile_Close(mfp);
  if (esl_msa_MinimGaps(msa, NULL, "-~") != eslOK) esl_fatal(msg);
  if (msa->alen        != 45)  esl_fatal(msg); /* orig =47, with one all - column and one all ~ column */
  if (msa->aseq[0][11] != 'L') esl_fatal(msg); /* L shifted from column 13->12 */
  if (msa->aseq[0][18] != 'T') esl_fatal(msg); /* T shifted from column 21->19 */
  esl_msa_Destroy(msa);

#ifdef eslAUGMENT_ALPHABET
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) esl_fatal(msg);
  if (esl_msafile_OpenDigital(abc, tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK) esl_fatal(msg);
  esl_msafile_Close(mfp);
  if (esl_msa_MinimGaps(msa, NULL, NULL) != eslOK) esl_fatal(msg);
  if (msa->alen        != 45)  esl_fatal(msg); /* orig =47, with one all - column and one all ~ column */
  if (esl_msa_Textize(msa) != eslOK) esl_fatal(msg);
  if (msa->aseq[0][11] != 'L') esl_fatal(msg); /* L shifted from column 13->12 */
  if (msa->aseq[0][18] != 'T') esl_fatal(msg); /* T shifted from column 21->19 */
  esl_msa_Destroy(msa);
  esl_alphabet_Destroy(abc);
#endif
  return;
}  

static void
utest_NoGaps(char *tmpfile)
{
  char        *msg = "NoGaps() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  ESL_MSA     *msa = NULL;
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET *abc = NULL;
#endif

  if (esl_msafile_Open(tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK)                                     esl_fatal(msg);
  esl_msafile_Close(mfp);
  if (esl_msa_NoGaps(msa, NULL, "-~") != eslOK) esl_fatal(msg);
  if (msa->alen        != 40)  esl_fatal(msg); /* orig =47, w/ 7 columns with gaps */
  if (msa->aseq[0][9]  != 'L') esl_fatal(msg); /* L shifted from column 13->10  */
  if (msa->aseq[0][16] != 'T') esl_fatal(msg); /* T shifted from column 21->17 */
  if (msa->aseq[0][39] != 'Y') esl_fatal(msg); /* Y shifted from column 47->40 */
  esl_msa_Destroy(msa);

#ifdef eslAUGMENT_ALPHABET
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) esl_fatal(msg);
  if (esl_msafile_OpenDigital(abc, tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK) esl_fatal(msg);
  esl_msafile_Close(mfp);
  if (esl_msa_NoGaps(msa, NULL, NULL) != eslOK) esl_fatal(msg);
  if (msa->alen        != 40)  esl_fatal(msg); /* orig =47, with one all - column and one all ~ column */
  if (esl_msa_Textize(msa) != eslOK) esl_fatal(msg);
  if (msa->aseq[0][9]  != 'L') esl_fatal(msg); /* L shifted from column 13->10  */
  if (msa->aseq[0][16] != 'T') esl_fatal(msg); /* T shifted from column 21->17 */
  if (msa->aseq[0][39] != 'Y') esl_fatal(msg); /* Y shifted from column 47->40 */
  esl_msa_Destroy(msa);
  esl_alphabet_Destroy(abc);
#endif
  return;
}  

static void
utest_SymConvert(char *tmpfile)
{
  char        *msg = "SymConvert() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  ESL_MSA     *msa = NULL;
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET *abc = NULL;
#endif

  if (esl_msafile_Open(tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK)                                     esl_fatal(msg);
  esl_msafile_Close(mfp);

  /* many->one version */
  if (esl_msa_SymConvert(msa, "VWY", "-")   != eslOK) esl_fatal(msg); /* 6 columns convert to all-gap: now 8/47 */
  if (esl_msa_MinimGaps(msa, NULL, "-~")    != eslOK) esl_fatal(msg); /* now we're 39 columns long */
  if (msa->alen                             != 39)    esl_fatal(msg);

  /* many->many version */
  if (esl_msa_SymConvert(msa, "DEF", "VWY") != eslOK) esl_fatal(msg);
  if (msa->aseq[0][4]                       != 'V')   esl_fatal(msg);
  if (msa->aseq[0][5]                       != 'W')   esl_fatal(msg);
  if (msa->aseq[0][23]                      != 'Y')   esl_fatal(msg); /* F in orig col 29; -5; converted to Y */

  /* bad calls */
#ifdef eslTEST_THROWING
  if (esl_msa_SymConvert(msa, "XXX", "XX")  != eslEINVAL) esl_fatal(msg); /* check for clean fail on mismatched args */
#endif
  esl_msa_Destroy(msa);
  
#ifdef eslAUGMENT_ALPHABET
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) esl_fatal(msg);
  if (esl_msafile_OpenDigital(abc, tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK) esl_fatal(msg);
  esl_msafile_Close(mfp);
#ifdef eslTEST_THROWING
  if (esl_msa_SymConvert(msa, "Tt", "Uu") != eslEINVAL) esl_fatal(msg); /* must cleanly fail on digital mode msa */
#endif
  esl_msa_Destroy(msa);
  esl_alphabet_Destroy(abc);
#endif
  return;
}

/* Exercise a boundary case: zero length MSA (alen=0) */
/* Given an input *digital* MSA as a starting point, we clone it, 
 * column subset it to zero length, then make sure that 
 * various MSA functions operate correctly on it;
 * then we textize it and test it in text mode; then we 
 * digitize it again, and throw it away.
 * (The input <msa> is unchanged.)
 */
static void
utest_ZeroLengthMSA(const char *tmpfile)
{
  char    *msg      = "zero length msa unit test failed";
  ESL_MSAFILE *mfp  = NULL;
  ESL_MSA *z1       = NULL;
  ESL_MSA *z2       = NULL;
  ESL_MSA *z3       = NULL;
  int     *useme    = NULL;
  int      nuseme   = 0;
  int      i;
  char     errbuf[eslERRBUFSIZE];


  /* Read a text mode alignment from the tmpfile */
  if (esl_msafile_Open(tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &z1) != eslOK)                                      esl_fatal(msg);
  esl_msafile_Close(mfp);

  /* make an alen=0 text alignment by column subsetting */
  nuseme = ESL_MAX(z1->alen, z1->nseq);
  if ((useme = malloc(sizeof(int) * nuseme)) == NULL)  esl_fatal(msg);
  for (i = 0; i < z1->alen; i++) useme[i] = 0;
  if (esl_msa_ColumnSubset(z1, errbuf, useme) != eslOK) esl_fatal(msg);

  /* These should all no-op if alen=0*/
  if (esl_msa_MinimGaps(z1, NULL, "-")!= eslOK) esl_fatal(msg);
  if (esl_msa_NoGaps(z1, NULL, "-")   != eslOK) esl_fatal(msg);
  if (esl_msa_SymConvert(z1,"RY","NN")!= eslOK) esl_fatal(msg);
  
  /* test sequence subsetting by removing the first sequence */
  for (i = 1; i < z1->nseq; i++) useme[i] = 1;  
  if (esl_msa_SequenceSubset(z1, useme, &z2) != eslOK) esl_fatal(msg);
  esl_msa_Destroy(z1);
  /* keep z2; we'll compare it to z3 in the end */
      
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET *abc;

  /* Now read the same alignment, in digital mode */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) esl_fatal(msg);
  if (esl_msafile_OpenDigital(abc, tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &z1) != eslOK) esl_fatal(msg);
  esl_msafile_Close(mfp);

  /* Now make an alen=0 alignment in digital mode */
  for (i = 0; i < z1->alen; i++) useme[i] = 0;
  if (esl_msa_ColumnSubset(z1, errbuf, useme) != eslOK) esl_fatal(msg);

  /* again these should all no-op if alen=0*/
  if (esl_msa_MinimGaps(z1, NULL, NULL) != eslOK) esl_fatal(msg);
  if (esl_msa_NoGaps(z1, NULL, NULL)    != eslOK) esl_fatal(msg);
  /* SymConvert throws EINVAL on a digital mode alignment */

  /* test sequence subsetting by removing the first sequence */
  for (i = 1; i < z1->nseq; i++) useme[i] = 1;  
  if (esl_msa_SequenceSubset(z1, useme, &z3) != eslOK) esl_fatal(msg);
  esl_msa_Destroy(z1);

  if ((z1 = esl_msa_Clone(z3))  == NULL)  esl_fatal(msg); /* z1 is now alen=0, digital */
  if (esl_msa_Textize(z3)       != eslOK) esl_fatal(msg); /* convert z3 back to text mode */
  if (esl_msa_Compare(z2, z3)   != eslOK) esl_fatal(msg); /* compare in text mode */
  if (esl_msa_Digitize(abc, z2) != eslOK) esl_fatal(msg); /* now z2 is digital */
  if (esl_msa_Compare(z1, z2)   != eslOK) esl_fatal(msg); /* compare digital mode z1,z2 */

  esl_alphabet_Destroy(abc);
  esl_msa_Destroy(z1);
  esl_msa_Destroy(z3);
#endif /*eslAUGMENT_ALPHABET*/

  esl_msa_Destroy(z2);
  free(useme);
}
#endif /* eslMSA_TESTDRIVE */
/*------------------------ end of unit tests --------------------------------*/



/*****************************************************************************
 * 10. Test driver
 *****************************************************************************/
#ifdef eslMSA_TESTDRIVE
/* 
 * gcc -g -Wall -o msa_utest -I. -DeslMSA_TESTDRIVE -DAUGMENT_KEYHASH esl_msa.c esl_keyhash.c easel.c -lm
 * gcc -g -Wall -o msa_utest -I. -DeslMSA_TESTDRIVE -DAUGMENT_ALPHABET esl_msa.c esl_alphabet.c easel.c -lm
 * gcc -g -Wall -o msa_utest -I. -DeslMSA_TESTDRIVE -DAUGMENT_SSI esl_msa.c esl_ssi.c easel.c -lm
 * gcc -g -Wall -o msa_utest -L. -I. -DeslMSA_TESTDRIVE esl_msa.c -leasel -lm
 * gcc -g -Wall -o msa_utest -L. -I. -DeslTEST_THROWING -DeslMSA_TESTDRIVE esl_msa.c -leasel -lm
 * ./msa_utest
 */
#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#ifdef eslAUGMENT_KEYHASH
#include "esl_keyhash.h"
#endif
#ifdef eslAUGMENT_RANDOM
#include "esl_random.h"
#endif
#ifdef eslAUGMENT_SSI
#include "esl_ssi.h"
#endif
#include "esl_msa.h"


int
main(int argc, char **argv)
{
  ESL_MSAFILE    *mfp  = NULL;
  ESL_MSA        *msa  = NULL;
  FILE           *fp   = NULL;
  char            tmpfile[16] = "esltmpXXXXXX"; /* tmpfile template */
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET   *abc  = NULL;
#endif

#ifdef eslTEST_THROWING
  esl_exception_SetHandler(&esl_nonfatal_handler);
#endif

  /* Create a known Stockholm test alignment in a tempfile.
   */
  if (esl_tmpfile_named(tmpfile, &fp) != eslOK) esl_fatal("failed to create tmpfile");
  write_known_msa(fp);
  fclose(fp);

  /* Read it back in for use in tests.
   */
  if (esl_msafile_Open(tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal("Failed to open MSA tmp file");
  if (esl_msa_Read(mfp, &msa)                                     != eslOK) esl_fatal("Failed to read MSA tmp file");
  esl_msafile_Close(mfp);

  /* Unit tests
   */
  utest_Create();
  utest_Destroy();
  utest_Expand();
  utest_Open(tmpfile);
  utest_Close(tmpfile);
  utest_Read(tmpfile);
  utest_Write(msa);
  utest_GuessFileFormat();
  utest_SequenceSubset(msa);
  utest_MinimGaps(tmpfile);
  utest_NoGaps(tmpfile);
  utest_SymConvert(tmpfile);
  utest_ZeroLengthMSA(tmpfile);	/* this tests in digital mode too if eslAUGMENT_ALPHABET */
  esl_msa_Destroy(msa);

#ifdef eslAUGMENT_ALPHABET
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) 
    esl_fatal("alphabet creation failed");
  if (esl_msafile_OpenDigital(abc, tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) 
    esl_fatal("MSA digital open failed");
  if (esl_msa_Read(mfp, &msa) != eslOK) 
    esl_fatal("MSA digital read failed");
  esl_msafile_Close(mfp);

  utest_CreateDigital(abc);
  utest_Digitize(abc, tmpfile);
  utest_Textize(abc, tmpfile);
  utest_OpenDigital(abc, tmpfile);
  utest_Write(msa);

  esl_alphabet_Destroy(abc);
  esl_msa_Destroy(msa);
#endif
  remove(tmpfile);
  exit(0);	/* success  */
}
#endif /*eslMSA_TESTDRIVE*/
/*-------------------- end of test driver ---------------------*/




/******************************************************************************
 * 11. Examples.
 *****************************************************************************/
/* The examples are also useful as i/o speed benchmarks. */



/* Example 1: 
 *   time ./example SSU_rRNA_5 > /dev/null
 *     without keyhash: [345.118u 31.564s 10:45.60 58.3%  SRE, Tue Sep  5 11:52:41 2006]
 *     with keyhash:    [33.353u 1.681s 0:35.04 99.9% SRE, Tue Sep  5 11:55:00 2006]
 */
#ifdef eslMSA_EXAMPLE
/*::cexcerpt::msa_example::begin::*/
/* An example of reading an MSA in text mode, and handling any returned errors.
   gcc -g -Wall -o example -I. -DeslMSA_EXAMPLE esl_msa.c easel.c 
   gcc -g -Wall -o example -I. -DeslMSA_EXAMPLE -DeslAUGMENT_KEYHASH esl_msa.c esl_keyhash.c easel.c
   ./example <MSA file>
 */
#include <stdio.h>
#include "easel.h"
#include "esl_msa.h"

int
main(int argc, char **argv)
{
  char        *filename = argv[1];
  int          fmt      = eslMSAFILE_UNKNOWN;
  ESL_MSAFILE *afp      = NULL;
  ESL_MSA     *msa      = NULL;
  int          nali     = 0;
  int          status;


  status = esl_msafile_Open(filename, fmt, NULL, &afp);
  if (status == eslENOTFOUND)    esl_fatal("Alignment file %s isn't readable\n", filename);
  else if (status == eslEFORMAT) esl_fatal("Couldn't determine format of %s\n",  filename);
  else if (status != eslOK)      esl_fatal("Alignment file open failed (error %d)\n", status);

  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;
      printf("alignment %5d: %15s: %6d seqs, %5d columns\n", 
	     nali, msa->name, msa->nseq, msa->alen);
      esl_msa_Write(stdout, msa, eslMSAFILE_STOCKHOLM);
      esl_msa_Destroy(msa);
    }
  if (status == eslEFORMAT)
    esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n",
	      afp->linenumber, afp->fname, afp->errbuf, afp->buf);
  else if (status != eslEOF)
    esl_fatal("Alignment file read failed with error code %d\n", status);

  esl_msafile_Close(afp);
  exit(0);
}
/*::cexcerpt::msa_example::end::*/
#endif /*eslMSA_EXAMPLE*/




#ifdef eslAUGMENT_ALPHABET
#ifdef eslMSA_EXAMPLE2
/*::cexcerpt::msa_example2::begin::*/
/* An example of reading an MSA in digital mode, after guessing
 * the alphabet by looking at the first alignment.
 *
   gcc -g -Wall -o example -I. -DeslMSA_EXAMPLE2 -DeslAUGMENT_ALPHABET esl_msa.c esl_alphabet.c easel.c 
   gcc -g -Wall -o example -I. -DeslMSA_EXAMPLE2 -DeslAUGMENT_ALPHABET -DeslAUGMENT_KEYHASH\
       esl_msa.c esl_keyhash.c esl_alphabet.c easel.c
   gcc -g -Wall -o example -I. -L. -DeslMSA_EXAMPLE2 esl_msa.c -leasel 
   ./example <MSA file>
 */
#include <stdio.h>
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"

int
main(int argc, char **argv)
{
  char         *filename = argv[1];
  int           fmt      = eslMSAFILE_UNKNOWN;
  int           type     = eslUNKNOWN;
  ESL_ALPHABET *abc      = NULL;
  ESL_MSAFILE  *afp      = NULL;
  ESL_MSA      *msa      = NULL;
  int           nali     = 0;
  int           status;

  /* First you open the msa file in normal text mode */
  status = esl_msafile_Open(filename, fmt, NULL, &afp);
  if (status == eslENOTFOUND)    esl_fatal("Alignment file %s isn't readable", filename);
  else if (status == eslEFORMAT) esl_fatal("Couldn't determine format of %s",  filename);
  else if (status != eslOK)      esl_fatal("Alignment file open failed (error code %d)", status);

  /* then you can guess the alphabet type - this looks at the first alignment */
  status = esl_msafile_GuessAlphabet(afp, &type);
  if      (status == eslEAMBIGUOUS) esl_fatal("Couldn't guess alphabet from first alignment in %s", filename);
  else if (status == eslEFORMAT)    esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nBad line is: %s\n",
					       afp->linenumber, afp->fname, afp->errbuf, afp->buf);
  else if (status == eslENODATA)    esl_fatal("Alignment file %s contains no data?", filename);
  else if (status != eslOK)         esl_fatal("Failed to guess alphabet (error code %d)\n", status);

  /* then you know how to create the alphabet */
  abc = esl_alphabet_Create(type);

  /* then you set the msafile into digital mode */
  esl_msafile_SetDigital(afp, abc);

  /* and now MSA's that you read are digital data in msa->ax, not text in msa->aseq */
  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;
      printf("alignment %5d: %15s: %6d seqs, %5d columns\n", 
	     nali, msa->name, msa->nseq, msa->alen);
      esl_msa_Write(stdout, msa, eslMSAFILE_STOCKHOLM);
      esl_msa_Destroy(msa);
    }
  if (status == eslEFORMAT)
    esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n",
	      afp->linenumber, afp->fname, afp->errbuf, afp->buf);
  else if (status != eslEOF)
    esl_fatal("Alignment file read failed with error code %d\n", status);

  esl_alphabet_Destroy(abc);
  esl_msafile_Close(afp);
  exit(0);
}
#endif /*eslMSA_EXAMPLE2*/
#endif /*eslAUGMENT_ALPHABET*/
/*::cexcerpt::msa_example2::end::*/
/*------------------------ end of examples -----------------------*/

 


/*****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/

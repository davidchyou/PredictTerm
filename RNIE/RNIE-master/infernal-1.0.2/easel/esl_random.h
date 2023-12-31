/* Easel's portable, threadsafe random number generator.
 * 
 * SRE, Wed Jul 14 11:23:57 2004 [St. Louis]
 * SVN $Id: esl_random.h 249 2008-04-24 19:19:50Z eddys $
 */
#ifndef ESL_RANDOM_INCLUDED
#define ESL_RANDOM_INCLUDED

typedef struct {
  long  seed;           /* reseed with this value, >0    */
  long  rnd1;           /* random number from LCG1       */
  long  rnd2;           /* random number from LCG2       */
  long  rnd;            /* random number we return       */
  long  tbl[64];        /* table for Bays/Durham shuffle */
  int   reseeding;	/* TRUE if seed is new           */
} ESL_RANDOMNESS;

/* esl_rnd_Roll(a) chooses a uniformly distributed integer
 * in the range 0..a-1, given an initialized ESL_RANDOMNESS r.
 */
#define esl_rnd_Roll(r, a)    ((int) (esl_random(r) * (a)))


/* 1. The ESL_RANDOMNESS object.
 */
extern ESL_RANDOMNESS *esl_randomness_Create(long seed);
extern ESL_RANDOMNESS *esl_randomness_CreateTimeseeded(void);
extern void            esl_randomness_Destroy(ESL_RANDOMNESS *r);
extern int             esl_randomness_Init(ESL_RANDOMNESS *r, long seed);
extern long            esl_randomness_GetSeed(const ESL_RANDOMNESS *r);

/* 2. The generator, esl_random().
 */
extern double esl_random(ESL_RANDOMNESS *r);

/* 3. Other fundamental sampling (including Gaussian, gamma).
 */
extern double esl_rnd_UniformPositive(ESL_RANDOMNESS *r);
extern double esl_rnd_Gaussian(ESL_RANDOMNESS *r, double mean, double stddev);
extern double esl_rnd_Gamma(ESL_RANDOMNESS *r, double a);

/* 4. Multinomial sampling from discrete probability n-vectors.
 */
extern int    esl_rnd_DChoose(ESL_RANDOMNESS *r, const double *p, int N);
extern int    esl_rnd_FChoose(ESL_RANDOMNESS *r, const float  *p, int N);


#endif /*ESL_RANDOM_INCLUDED*/

/*****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/

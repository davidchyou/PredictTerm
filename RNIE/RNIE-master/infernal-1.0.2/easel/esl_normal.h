/* Statistical routines for normal distributions
 * 
 * SRE, Tue Nov 21 14:29:02 2006 [Janelia]
 * SVN $Id: esl_normal.h 269 2008-06-19 13:47:41Z eddys $
 */

#ifndef ESL_NORMAL_INCLUDED
#define ESL_NORMAL_INCLUDED

extern double esl_normal_pdf   (double x, double mu, double sigma);
extern double esl_normal_logpdf(double x, double mu, double sigma);
extern double esl_normal_cdf   (double x, double mu, double sigma);
extern double esl_normal_surv  (double x, double mu, double sigma);

#endif /*ESL_NORMAL_INCLUDED*/
/*****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/

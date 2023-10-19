/* esl_msaweight.h
 * Sequence weighting algorithms.
 * 
 * SVN $Id: esl_msaweight.h 344 2009-06-15 15:32:24Z nawrockie $
 * SRE, Sun Nov  5 09:11:13 2006 [Janelia]
 */
#ifndef ESL_MSAWEIGHT_INCLUDED
#define ESL_MSAWEIGHT_INCLUDED

#include <esl_msa.h>

extern int esl_msaweight_GSC(ESL_MSA *msa);
extern int esl_msaweight_PB(ESL_MSA *msa);
extern int esl_msaweight_BLOSUM(ESL_MSA *msa, double maxid);
extern int esl_msaweight_IDFilter(const ESL_MSA *msa, double maxid, ESL_MSA **ret_newmsa);


#endif /*ESL_MSAWEIGHT_INCLUDED*/

/*****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/

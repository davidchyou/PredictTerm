/* scancyk.c
 * SRE, Thu May  2 11:50:48 2002 [AA 3050 SFO->STL]
 * SVN $Id: scancyk.c 2244 2007-12-05 13:22:15Z nawrockie $
 * 
 * CYK alignment: multihit, local, database scanning mode.
 * [xref STL6 p47]
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

#include "easel.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

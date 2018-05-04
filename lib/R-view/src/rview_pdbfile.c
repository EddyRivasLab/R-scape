/* rview_pddfile - functions to read a pdbx/mmcif file
 * Contents:
 *
 * ER, Thu May  3 19:21:13 EDT 2018 [Verril Farm] 
 * SVN $Id:$
 */

#include "rview_config.h"

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include <math.h>
#include <float.h>
	
#include "easel.h"

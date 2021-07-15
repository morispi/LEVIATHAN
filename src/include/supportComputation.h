#ifndef __MYVC_SUPPORT_COMPUTATION__
#define __MYVC_SUPPORT_COMPUTATION__

#include "SVSupports.h"
#include <vector>
#include "api/BamAlignment.h"

using namespace std;
using namespace BamTools;

/**
	Compute the SV support of a given candidate.
	@param commonAlignments1 the alignments of the 1st region that have a matching read in the 2nd region
	@param commonAlignments2 the alignments of the 2nd region that have a matching read in the 1st region
	@param begR1 the beginning position of region1
	@param begR2 the beginning position of region2
	@param posSupportRegion1 an array representing the split reads positions support of the first region
	@param posSupportRegion2 an array representing the split reads positions support of the second region
	@return a structure containing all support information
*/
SVSupports computeSVSupport(vector<BamAlignment>& commonAlignments1, vector<BamAlignment>& commonAlignments2, unsigned begR1, unsigned begR2, unsigned posSupportRegion1[], unsigned posSupportRegion2[]);

#endif
#ifndef __MYVC_SUPPORT_COMPUTATION__
#define __MYVC_SUPPORT_COMPUTATION__

#include "SVSupports.h"
#include <vector>
#include "api/BamAlignment.h"

using namespace std;
using namespace BamTools;

SVSupports computeSVSupport(vector<BamAlignment>& commonAlignments1, vector<BamAlignment>& commonAlignments2, unsigned begR1, unsigned begR2, unsigned posSupportRegion1[], unsigned posSupportRegion2[]);

#endif
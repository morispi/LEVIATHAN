#ifndef __MYVC_SV_VALIDATION__
#define __MYVC_SV_VALIDATION__

#include "robin_hood.h"
#include <string>
#include "misc.h"
#include "SVSupports.h"
#include "StructuralVariant.h"

/**
	Validate the SV according to their support information
	@param calledSVs a map associating region pairs to their SV support, for SVs detected previously
	@return a set of valid SVs
*/
robin_hood::unordered_set<StructuralVariant> validatesSVs(robin_hood::unordered_map<pair<string*, string*>, SVSupports, hashPointersPairs>& calledSVs);

#endif
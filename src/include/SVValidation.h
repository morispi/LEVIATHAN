#ifndef __MYVC_SV_VALIDATION__
#define __MYVC_SV_VALIDATION__

#include "robin_hood.h"
#include <string>
#include "misc.h"
#include "SVSupports.h"
#include "StructuralVariant.h"

robin_hood::unordered_set<StructuralVariant> validatesSVs(robin_hood::unordered_map<pair<string*, string*>, SVSupports, hashPointersPairs>& calledSVs);

#endif
#ifndef __MYVC_SV_PROCESSING__
#define __MYVC_SV_PROCESSING__

#include "robin_hood.h"
#include <string>
#include "misc.h"
#include "SVSupports.h"

robin_hood::unordered_set<string> validatesSVs(robin_hood::unordered_map<pair<string*, string*>, SVSupports, hashPointersPairs>& calledSVs);

void outputSVs(robin_hood::unordered_set<string>& finalSVs);

#endif
#ifndef __MYVC_SV_OUTPUT__
#define __MYVC_SV_OUTPUT__

#include "robin_hood.h"
#include <string>
#include "misc.h"
#include "SVSupports.h"
#include "StructuralVariant.h"

void outputSVsAsBed(robin_hood::unordered_set<StructuralVariant>& finalSVs);

void outputSVsAsVCF(robin_hood::unordered_set<StructuralVariant>& finalSVs, string cmdLine, string reference);

#endif
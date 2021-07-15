#ifndef __MYVC_SV_OUTPUT__
#define __MYVC_SV_OUTPUT__

#include "robin_hood.h"
#include <string>
#include "misc.h"
#include "SVSupports.h"
#include "StructuralVariant.h"

/**
	Output SVs in bed format, in stdout
	@param finalSVs the set of SVs to output
*/
void outputSVsAsBed(robin_hood::unordered_set<StructuralVariant>& finalSVs);

/**
	Output SVs in VCF format
	@param finalSVs the set of SVs to output
	@param cmdLine the command line that was used to run LEVIATHAN
	@param reference the name of the reference genome file
	@param outputFile file where to output the SVs
*/
void outputSVsAsVCF(robin_hood::unordered_set<StructuralVariant>& finalSVs, string cmdLine, string reference, string outputFile);

#endif
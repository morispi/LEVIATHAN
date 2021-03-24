#ifndef __MYVC_ALIGNMENTS_PROCESSING__
#define __MYVC_ALIGNMENTS_PROCESSING__

#include <string>
#include <vector>
#include "api/BamAlignment.h"
#include "robin_hood.h"
#include "SVSupports.h"

using namespace std;
using namespace BamTools;

pair<robin_hood::unordered_set<string>, vector<BamAlignment>> extractAlignmentsAndHeadersFromRegion(string& bamFile, string& region, unsigned additionalSize);

vector<BamAlignment> retrieveAlignmentsWithCommonHeadersAndTreatInsertions(vector<BamAlignment>& alignments, robin_hood::unordered_set<string>& headers, SVSupports& support, unsigned posSupportRegion[], int32_t beg);

pair<robin_hood::unordered_set<string>, vector<BamAlignment>> extractAlignmentsWithCommonHeadersAndTreatInsertions(string& bamFile, string& region, unsigned additionalSize, robin_hood::unordered_set<string>& headers, SVSupports& support, unsigned posSupportRegion[], int32_t beg);

#endif
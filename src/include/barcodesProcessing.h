#ifndef __MYVC_BARCODES_PROCESSING__
#define __MYVC_BARCODES_PROCESSING__

#include <string>
#include <vector>
#include "robin_hood.h"
#include "indexManagementBam.h"
#include "misc.h"

using namespace std;

robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> processBarcode(int id, pair<int32_t, int32_t> beg, pair<int32_t, int32_t> end, string& bamFile, BarcodesPositionsIndex& barcodesPositionsIndex, const barcode& b, int barcodesSize, int& distance);

robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> processBarcodes(int nbThreads, robin_hood::unordered_map<string, int>& refIDs, vector<string>& regionsList, unsigned nbBins, string& bamFile, BarcodesPositionsIndex& barcodesPositionsIndex, int barcodesSize, int distance);

Thresholds analyzeDistribution(robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates);

#endif

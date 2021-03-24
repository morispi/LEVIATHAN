#ifndef __MYVC_GLOBAL_VARS__
#define __MYVC_GLOBAL_VARS__

#include "robin_hood.h"
#include <string>
#include <mutex>
#include <vector>

using namespace std;

// Global variables used throughout the functions
extern mutex countMtx;
extern unsigned processedBarcodes;
extern unsigned totalBarcodes;
extern unsigned processedCandidates;
extern unsigned totalCandidates;
extern robin_hood::unordered_map<int, vector<string>> windows;
extern vector<string> regionsList;
extern robin_hood::unordered_map<string, int32_t> refIDs;
extern unsigned windowSize;
extern robin_hood::unordered_map<std::string, std::vector<bool>> genomeIndex;
extern unsigned mediumVariantsSize;
extern unsigned largeVariantsSize;
extern double smallVariantsRate;
extern double mediumVariantsRate;
extern double largeVariantsRate;
extern unsigned duplicatesDistance;
extern unsigned poolSize;
extern unsigned minBarcodes;

#endif
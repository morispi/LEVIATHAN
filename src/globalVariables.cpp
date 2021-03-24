#include "globalVariables.h"

mutex countMtx;
unsigned processedBarcodes;
unsigned totalBarcodes;
unsigned processedCandidates;
unsigned totalCandidates;
robin_hood::unordered_map<int, vector<string>> windows;
vector<string> regionsList;
robin_hood::unordered_map<string, int32_t> refIDs;
unsigned windowSize;
robin_hood::unordered_map<std::string, std::vector<bool>> genomeIndex;
unsigned mediumVariantsSize;
unsigned largeVariantsSize;
double smallVariantsRate;
double mediumVariantsRate;
double largeVariantsRate;
unsigned duplicatesDistance;
unsigned poolSize;
unsigned minBarcodes;
#ifndef __MYVC_BARCODES_PROCESSING__
#define __MYVC_BARCODES_PROCESSING__

#include <string>
#include <vector>
#include "robin_hood.h"
#include "indexManagementBam.h"
#include "misc.h"

using namespace std;

/**
	Process a given barcode to retrieve region pairs that share it.
	@param id identifier for running threads
	@param beg start of the index interval to consider
	@param end end of the index interval to consider
	@param bamFile BAM file containing the alignments
	@param barcodesPositionsIndex barcodes index built using LRez
	@param barcode to process
	@param barcodesSize size of the index
	@param minVariantSize minimum variant size to consider. Region pairs closer than this value will not be considered
	@param skipTranslocations whether to process translocation SVs or not
	@return a map associating a region to the list of regions it shares the barcode with, along with the number of barcodes shared by the regions
*/
robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> processBarcode(int id, pair<int32_t, int32_t> beg, pair<int32_t, int32_t> end, string& bamFile, BarcodesPositionsIndex& barcodesPositionsIndex, const barcode& b, int barcodesSize, int& distance, bool skipTranslocations);

/**
	Process all the barcodes and compute the number of common barcodes between all region pairs
	@param nbThreads number of threads to use
	@param refIDs map associating a chromosome name to its integer value
	@param regionsList list containing all regions of the reference genome
	@param nbBins number of bins to separate the index into
	@param bamFile BAM file containing the alignments
	@param barcodesPositionsIndex barcodes index built using LRez
	@param barcodesSize size of the index
	@param minVariantSize minimum variant size to consider. Region pairs closer than this value will not be considered
	@param skipTranslocations whether to process translocation SVs or not
	@return a map associating a region to the list of regions it shares barcodes with, along with the number of barcodes shared by the regions
*/
robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> processBarcodes(int nbThreads, robin_hood::unordered_map<string, int>& refIDs, vector<string>& regionsList, unsigned nbBins, string& bamFile, BarcodesPositionsIndex& barcodesPositionsIndex, int barcodesSize, int distance, bool skipTranslocations);

Thresholds analyzeDistribution(robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates);

#endif

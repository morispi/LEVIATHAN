#ifndef __MYVC_ALIGNMENTS_PROCESSING__
#define __MYVC_ALIGNMENTS_PROCESSING__

#include <string>
#include <vector>
#include "api/BamAlignment.h"
#include "robin_hood.h"
#include "SVSupports.h"

using namespace std;
using namespace BamTools;

/**
	Extract alignments and headers of a BAM file given a genomic region.
	@param bamFile BAM file to extract from
	@param region region of interest
	@param additionalSize extension size on the left and on the right of the region of interest
	@return a pair containing the set of headers and the vector of alignments
*/
pair<robin_hood::unordered_set<string>, vector<BamAlignment>> extractAlignmentsAndHeadersFromRegion(string& bamFile, string& region, unsigned additionalSize);

/**
	Retrieve alignments having their headers in a specified list, and treat insertion variants.
	@param alignments vector of alignments to process
	@param headers list of headers of interest
	@param support structure containing support information, filled with insertion supports
	@param posSupportRegion[] an array filled with the split read support information for insertions
	@param beg beginning position of the processed region
	@return a vector containing alignments having their headers in the specified set
*/
vector<BamAlignment> retrieveAlignmentsWithCommonHeadersAndTreatInsertions(vector<BamAlignment>& alignments, robin_hood::unordered_set<string>& headers, SVSupports& support, unsigned posSupportRegion[], int32_t beg);

/**
	Extract alignments and headers of a BAM file given a genomic region, if their headers are in a specified list and treat insertion variants.
	@param bamFile BAM file to extract from
	@param region region of interest
	@param additionalSize extension size on the left and on the right of the region of interest
	@param headers list of headers of interest
	@param support structure containing support information, filled with insertion supports
	@param posSupportRegion[] an array filled with the split read support information for insertions
	@param beg beginning position of the processed region
	@return a pair containing the set of headers and the vector of alignments having their headers in the specified set
*/
pair<robin_hood::unordered_set<string>, vector<BamAlignment>> extractAlignmentsWithCommonHeadersAndTreatInsertions(string& bamFile, string& region, unsigned additionalSize, robin_hood::unordered_set<string>& headers, SVSupports& support, unsigned posSupportRegion[], int32_t beg);

#endif
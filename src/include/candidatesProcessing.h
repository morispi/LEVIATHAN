#ifndef __MYVC_CANDIDATES_PROCESSING__
#define __MYVC_CANDIDATES_PROCESSING__

#include <vector>
#include <string>
#include "SVSupports.h"
#include "robin_hood.h"
#include "misc.h"

using namespace std;

/**
	Process a set SV candidates, defined by a region and the list of regions it shares barcodes with.
	@param id identifier for running threads
	@param candidate candidate of interest, pair associating the region of interest to the list of regions it shares barcodes with, along with their number of shared barcodes
	@bamFile BAM file containing the alignments
	@return a vector associating pairs of regions to the SV support they share
*/
vector<pair<pair<string*, string*>, SVSupports>> processCandidate(int id, pair<string*, robin_hood::unordered_map<string*, unsigned>>& candidate, string& bamFile);

/**
	Process all SV candidates, defined by a region and the list of regions it shares barcodes with.
	@param nbThreads number of threads to use
	@param candidates map containing all SV candidates, associating a given region to the list of regions it shares barcodes with, along with their number of shared barcodes
	@bamFile BAM file containing the alignments
	@return a map associating pairs of regions to the SV support they share
*/
robin_hood::unordered_map<pair<string*, string*>, SVSupports, hashPointersPairs> processCandidates(unsigned nbThreads, robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates, string& bamFile);

/**
	Remove candidates that share an insufficient number of barcodes
	@param candidates map containing all SV candidates, associating a given region to the list of regions it shares barcodes with, along with their number of shared barcodes
	@param th structure containing the thresholds values for close / moderately distant / distant region pairs
	@param maxRegionsLinks maximum number of regions a given region can be associated to, candidates associated to more regions will be filtered out

*/
void removeInvalidCandidates(robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates, Thresholds th, unsigned maxRegionsLinks);

/**
	Remove candidates that are associated to too many regions
	@param candidates map containing all SV candidates, associating a given region to the list of regions it shares barcodes with, along with their number of shared barcodes
	@param maxRegionsLinks maximum number of regions a given region can be associated to, candidates associated to more regions will be filtered out

*/
void removeCandidatesRegionsLinks(robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates, unsigned maxRegionsLinks);

/**
	Save all SV candidates to a file
	@param file file where to save SV candidates
	@param candidates map containing all SV candidates, associating a given region to the list of regions it shares barcodes with, along with their number of shared barcodes

*/
void saveSVCandidates(string file, robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> candidates);

/**
	Load SV candidates from a file
	@param file file containing SV candidates
	@return a map containing all SV candidates, associating a given region to the list of regions it shares barcodes with, along with their number of shared barcodes

*/
robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> loadSVCandidates(string file);

#endif
#ifndef __MYVC_MISC__
#define __MYVC_MISC__

#include <string>
#include <vector>
#include "api/BamAlignment.h"
#include "indexManagementBam.h"
#include "SVSupports.h"

/**
	Structure defining the thresholds for close / moderately distant and distant region pairs.
*/
struct Thresholds {
    unsigned closeTh;
    unsigned averageTh;
    unsigned farTh;
};

/**
	Structure defining a function allowing to hash pairs of pointers.
*/
struct hashPointersPairs { 
    template <class T1, class T2> 
    size_t operator()(const pair<T1*, T2*>& p) const
    { 
        auto hash1 = hash<T1>{}(*(p.first)); 
        auto hash2 = hash<T2>{}(*(p.second)); 
        return hash1 ^ hash2; 
    } 
};

/**
	Structure defining a function allowing to hash pairs.
*/
struct hashPairs { 
    template <class T1, class T2> 
    size_t operator()(const pair<T1, T2>& p) const
    { 
        auto hash1 = hash<T1>{}(p.first); 
        auto hash2 = hash<T2>{}(p.second); 
        return hash1 ^ hash2; 
    } 
};

/**
	Checks wether a file is empty or not.
	@param fileName file to check
	@return true if the file is empty, false otherwise
*/
bool isEmptyFile(string fileName);

/**
	Split a string according to a delimiter.
	@param s string to split
	@param delimiter to use
	@return a vector containing the different splits
*/
vector<string> splitString(string s, string delimiter);

/**
	Return the chromosome and beginning position of a region.
	@param region region of interest in format chr:beg-end
	@return a pair containing the name of the chromosome and the begining position of the region
*/
pair<string, int32_t> regionChrAndBegPosition(string region);

/**
	Return the beginning position of a region.
	@param region region of interest in format chr:beg-end
	@return the begining position of the region
*/
int32_t regionBegPosition(string region);


/**
	Extract the regions of a given chromosome.
	@param contig chromosome of interest
	@param id identifier of the chromosome
	@param contigSize size of the chromosome
*/
vector<string> extractWindowsRegions(string contig, unsigned id, int32_t contigSize);

/**
	Return the region corresponding to a given position on a chromosome.
	@param regions list of regions of the chromosome of interest
	@param position of interest
	@return the region corresponding to the position, in format chr:beg-end

*/
string* positionToRegion(vector<string>& regions, int32_t position);

/**
	Compare two alignments according to the name of the chromosome
	@param al1 first alignment to compare
	@param al2 second alignment to compare
	@return a positive value is al1 < al2 and a negative value otherwise
*/
bool orderBamAlignments (BamAlignment al1, BamAlignment al2);

/**
	Compute support of the strongest SV evidence contained in the structure.
	@param s the structure containing all support information
	@return the support of the strongest SV evidence
*/
unsigned computeMaxSupport(SVSupports s);

/**
	Compute the support of discordant paired reads contained in the structure.
	@param s the structure containing all support information
	@return the support of discordant paired reads
*/
unsigned computePaiedReadsSupport(SVSupports s);

/**
	Compare the SV support of two pairs of regions.
	@param a the first region pair / support to compare
	@param a the second region pair / support to compare
	@return a positive value is the support of a is greater than that of b, and a negative value otherwise
*/
bool compareSVSupport(const pair<pair<string*, string*>, SVSupports>& a, const pair<pair<string*, string*>, SVSupports>& b);

/**
	Remove the redundancy of the index, ie sequences of occurrence position that correspond to the same region.
	In that case, only one occurrence position is kept.
	@param barcodesPositionsIndex the index to remove redundancy from
*/
void removeIndexRedundancy(BarcodesPositionsIndex& barcodesPositionsIndex);

/**
	Prepare the auxiliary data necessary to LEVIATHAN, ie the list of all possible regions, the map associating a chromosome name to its ID, etc
	@param bamFile the BAM file to prepare data for
*/
void prepareAuxiliaryData(string& bamFile);

/**
	Translate a string into a BamRegion
	@param reader a BamReader open on the desired BAM file
	@param s the string to translate
	@return a BamRegion corresponding to the input string
*/
BamRegion stringToBamRegion(BamReader& reader, string s);

#endif
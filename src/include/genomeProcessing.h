#ifndef __MYVC_GENOME_PROCESSING__
#define __MYVC_GENOME_PROCESSING__

#include "robin_hood.h"
#include <string>
#include <vector>

using namespace std;

string twoBitsToString(bool b1, bool b2);

/**
	Convert a string into 2 bits per nucleotide representation.
	@param str string to covert
	@return the string in 2 bits per nucleotide representation
*/
vector<bool> fullstr2num(const string& str);

/**
	Convert a string in 2 bits per nucleotide representation into its string representation.
	@param num the string in 2 bits per nucleotide representation
	@return the string string representation
*/
string fullnum2str(vector<bool> num);

/**
	Index the reference genome, with strings in 2 bits per nuclotide representation.
	@param file the file containing the reference genome
	@return a map associating a chromosome name to it's sequence in 2 bits per nucleotide representation
*/
robin_hood::unordered_map<string, vector<bool>> indexGenome(string file);

/**
	Extract a given region of the reference genome
	@param region the region of interest in format chr:beg-end
	@return the sequence corresping to the region of interest
*/
string extractGenomeRegion(string region);

#endif
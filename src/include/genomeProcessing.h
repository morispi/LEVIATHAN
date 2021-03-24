#ifndef __MYVC_GENOME_PROCESSING__
#define __MYVC_GENOME_PROCESSING__

#include "robin_hood.h"
#include <string>
#include <vector>

using namespace std;

vector<bool> fullstr2num(const string& str);

string fullnum2str(vector<bool> num);

robin_hood::unordered_map<string, vector<bool>> indexGenome(string readsFile);

string extractGenomeRegion(string region);

#endif
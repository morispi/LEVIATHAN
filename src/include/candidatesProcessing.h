#ifndef __MYVC_CANDIDATES_PROCESSING__
#define __MYVC_CANDIDATES_PROCESSING__

#include <vector>
#include <string>
#include "SVSupports.h"
#include "robin_hood.h"
#include "misc.h"

using namespace std;

vector<pair<pair<string*, string*>, SVSupports>> processCandidate(int id, pair<string*, robin_hood::unordered_map<string*, unsigned>>& candidate, string& bamFile);

robin_hood::unordered_map<pair<string*, string*>, SVSupports, hashPointersPairs> processCandidates(unsigned nbThreads, robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates, string& bamFile);

void removeInvalidCandidates(robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates, Thresholds th, unsigned maxRegionsLinks);

void removeCandidatesRegionsLinks(robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates, unsigned maxRegionsLinks);

void saveSVCandidates(string file, robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> candidates);

robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> loadSVCandidates(string file);

#endif
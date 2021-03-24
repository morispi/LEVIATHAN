#ifndef __MYVC_MISC__
#define __MYVC_MISC__

#include <string>
#include <vector>
#include "api/BamAlignment.h"
#include "indexManagementBam.h"
#include "SVSupports.h"

struct Thresholds {
    unsigned closeTh;
    unsigned averageTh;
    unsigned farTh;
};

struct hashPointersPairs { 
    template <class T1, class T2> 
    size_t operator()(const pair<T1*, T2*>& p) const
    { 
        auto hash1 = hash<T1>{}(*(p.first)); 
        auto hash2 = hash<T2>{}(*(p.second)); 
        return hash1 ^ hash2; 
    } 
};

struct hashPairs { 
    template <class T1, class T2> 
    size_t operator()(const pair<T1, T2>& p) const
    { 
        auto hash1 = hash<T1>{}(p.first); 
        auto hash2 = hash<T2>{}(p.second); 
        return hash1 ^ hash2; 
    } 
};

bool isEmptyFile(string fileName);

vector<string> splitString(string s, string delimiter);

pair<string, int32_t> regionChrAndBegPosition(string region);

int32_t regionBegPosition(string region);

vector<string> extractWindowsRegions(string contig, unsigned id, int32_t contigSize);

string* positionToRegion(vector<string>& regions, int32_t position);

bool orderBamAlignments (BamAlignment al1, BamAlignment al2);

unsigned computeMaxSupport(SVSupports s);

unsigned computePaiedReadsSupport(SVSupports s);

bool compareSVSupport(const pair<pair<string*, string*>, SVSupports>& a, const pair<pair<string*, string*>, SVSupports>& b);

void removeIndexRedundancy(BarcodesPositionsIndex& barcodesPositionsIndex);

void prepareAuxiliaryData(string& bamFile);

BamRegion stringToBamRegion(BamReader& reader, string s);

#endif
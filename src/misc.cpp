#include "misc.h"
#include "globalVariables.h"

bool isEmptyFile(string fileName) {
    std::ifstream file(fileName);
    bool res = file.peek() == std::ifstream::traits_type::eof();
    file.close();
    return res;
}

vector<string> splitString(string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

int32_t regionBegPosition(string region) {
    vector<string> t = splitString(region, ":");
    vector<string> tt = splitString(t[1], "-");

    return static_cast<uint32_t>(stoul(tt[0]));
}

pair<string, int32_t> regionChrAndBegPosition(string region) {
    vector<string> t = splitString(region, ":");
    vector<string> tt = splitString(t[1], "-");

    return make_pair(t[0], static_cast<uint32_t>(stoul(tt[0])));
}

vector<string> extractWindowsRegions(string contig, unsigned id, int32_t contigSize) {
	vector<string> res;
	for (int32_t i = 0; i < contigSize; i += windowSize) {
		res.push_back(contig + ":" + to_string(i) + "-" + to_string(i + windowSize - 1));
	}
	res.push_back(contig + ":" + to_string(contigSize - windowSize + 1) + "-" + to_string(contigSize));

	return res;
}

string* positionToRegion(vector<string>& regions, int32_t position) {
	return &regions[position / windowSize];
}

bool orderBamAlignments (BamAlignment al1, BamAlignment al2) { 
    return al1.Name < al2.Name;
}

unsigned computeMaxSupport(SVSupports s) {
    return max(s.deletion, max(s.duplication, max(s.inversion, max(s.insertion, s.translocation))));
}

unsigned computePairedReadsSupport(SVSupports s) {
    return s.deletion + s.duplication + s.inversion + s.insertion + s.translocation;
}

bool compareSVSupport(const pair<pair<string*, string*>, SVSupports>& a, const pair<pair<string*, string*>, SVSupports>& b) { 
    // return a.second.barcodes + computeMaxSupport(a.second) > b.second.barcodes + computeMaxSupport(b.second); 
    // return a.second.barcodes > b.second.barcodes;
    return a.second.barcodes + computePairedReadsSupport(a.second) > b.second.barcodes + computePairedReadsSupport(b.second); 
}

void removeIndexRedundancy(BarcodesPositionsIndex& barcodesPositionsIndex) {
    vector<barcode> vbarcodes(barcodesPositionsIndex.size());
    unsigned i = 0;
    for (auto p : barcodesPositionsIndex) {
        vbarcodes[i] = p.first;
        i++;
    }
    for (auto b : vbarcodes) {
        vector<pair<int32_t, int32_t>> v;
        v.push_back(barcodesPositionsIndex[b][0]);
        for (unsigned i = 1; i < barcodesPositionsIndex[b].size(); i++) {
            if (positionToRegion(windows[barcodesPositionsIndex[b][i].first], barcodesPositionsIndex[b][i].second) != positionToRegion(windows[barcodesPositionsIndex[b][i - 1].first], barcodesPositionsIndex[b][i - 1].second)) {
                v.push_back(barcodesPositionsIndex[b][i]);
            }
        }
        barcodesPositionsIndex[b] = v;
    }
}

void prepareAuxiliaryData(string& bamFile) {
    // Open BAM file and check if it has an existing index
    BamReader reader;
    if (!reader.Open(bamFile)) {
        fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
        exit(EXIT_FAILURE);
    }
    if (!reader.LocateIndex()) {
        fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
        exit(EXIT_FAILURE);
    }

    // Get a vector containing reference sequences data
    RefVector rv = reader.GetReferenceData();

    for (RefData d : rv) {
        int id = reader.GetReferenceID(d.RefName);
        BamAlignment al;
        if (id == -1) {
            fprintf(stderr, "Cannot find refence with ID %s.\n", d.RefName.c_str());
            exit(EXIT_FAILURE);
        }   
        // Only process the chromosome if it has alignments
        if (!reader.SetRegion(id, 0, id, d.RefLength - 1)) {
            fprintf(stderr, "Error while attempting to jump to region.\n");
            exit(EXIT_FAILURE);
        }
        if (reader.GetNextAlignment(al)) {
            vector<string> w = extractWindowsRegions(d.RefName, id, d.RefLength);
            for (string ww : w) {
                regionsList.push_back(ww);
            }
            windows[reader.GetReferenceID(d.RefName)] = w;
            refIDs[d.RefName] = reader.GetReferenceID(d.RefName);
        }
    }
    reader.Close();
}

BamRegion stringToBamRegion(BamReader& reader, string s) {
        BamRegion r;

        vector<string> t = splitString(s, ":");
    if (t.size() != 2) {
                return r;
        }

    vector<string> p = splitString(t[1], "-");
    if (p.size() != 2) {
        return r;
    }

        int leftID = reader.GetReferenceID(t[0]);
        if (leftID == -1) {
                fprintf(stderr, "Cannot find refence with ID %s.\n", t[0].c_str());
                exit(EXIT_FAILURE);
        }

        BamRegion res(leftID, stoi(p[0]), leftID, stoi(p[1]));
        return res;
}
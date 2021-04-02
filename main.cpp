#include <set>
#include "string.h"
#include <vector>
#include <unordered_map>
#include "barcodesExtraction.h"
#include "barcodesComparison.h"
#include "alignmentsRetrieval.h"
#include <mutex>
#include <future>
#include "CTPL/ctpl_stl.h"

// Global variables used throughout the functions
mutex countMtx;
unsigned processedBarcodes;
unsigned totalBarcodes;
unsigned processedCandidates;
unsigned totalCandidates;
robin_hood::unordered_map<int, vector<string>> windows;
vector<string> regionsList;
robin_hood::unordered_map<string, int32_t> refIDs;
unsigned windowSize;

struct SVSupports {
	unsigned barcodes = 0;
	unsigned alignments1 = 0;
	unsigned alignments2 = 0;
	unsigned deletion = 0;
	unsigned duplication = 0;
	unsigned inversion = 0;
	unsigned insertion = 0;
	unsigned translocation = 0;
	unsigned splitReads = 0;
	unsigned breakpoint1 = 0;
	unsigned support1 = 0;
	unsigned breakpoint2 = 0;
	unsigned support2 = 0;
};

void saveSVCandidates(string file, robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> candidates) {
	ofstream out;
	out.open(file);
	if (!out.is_open()) {
		fprintf(stderr, "Unable to open file %s.", file.c_str());
		exit(EXIT_FAILURE);
	}

	for (auto c : candidates) {
		vector<string> r1 = splitString(*(c.first), ":");
		vector<string> rr1 = splitString(r1[1], "-");
		for (auto cc : c.second) {
			vector<string> r2 = splitString(*(cc.first), ":");
			vector<string> rr2 = splitString(r2[1], "-");

			out << r1[0] << "\t" << rr1[0] << "\t" << rr1[1] << "\t" << r2[0] << "\t" << rr2[0] << "\t" << rr2[1] << "\t" << cc.second << endl;
		}
	}

	out.close();
}

int32_t regionBegPosition(string region) {
	vector<string> t = splitString(region, ":");
    vector<string> tt = splitString(t[1], "-");

    return static_cast<uint32_t>(stoul(tt[0]));
}

struct hashPointersPairs { 
    template <class T1, class T2> 
    size_t operator()(const pair<T1*, T2*>& p) const
    { 
        auto hash1 = hash<T1>{}(*(p.first)); 
        auto hash2 = hash<T2>{}(*(p.second)); 
        return hash1 ^ hash2; 
    } 
};

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


// TODO TOUT PASSER EN REF AUX THREADS
robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> processBarcode(int id, pair<int32_t, int32_t> beg, pair<int32_t, int32_t> end, string& bamFile, BarcodesPositionsIndex& barcodesPositionsIndex, barcode& b, int barcodesSize, int& distance) {
	robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> candidates;
	auto v = barcodesPositionsIndex[b];

	// Search for pairs of sufficiently distant windows (reg1, reg2) that share common barcodes
	unsigned i = 0;
	while (i < v.size() and (v[i].first < beg.first or (v[i].first == beg.first and v[i].second < beg.second))) {
		i++;
	}
	while (i < v.size() and (v[i].first < end.first or (v[i].first == end.first and v[i].second < end.second))) { 
		for (unsigned j = i + 1; j < v.size(); j++) {
			string* reg1 = positionToRegion(windows[v[i].first], v[i].second);
			string* reg2 = positionToRegion(windows[v[j].first], v[j].second);
			unsigned begR1 = (v[i].second / windowSize) * windowSize;
			unsigned begR2 = (v[j].second / windowSize) * windowSize;
			if (v[i].first != v[j].first) {
				if (v[i].first < v[j].first) {
					candidates[reg1][reg2]++;
				} else {
					candidates[reg2][reg1]++;
				}
			// TODO: Vérifier la distance entre les fenêtres
			// distance was 2 * windowSize
			} else if (reg1 != reg2 and begR2 >= begR1 + distance) {
				candidates[reg1][reg2]++;
			} else if (reg1 != reg2 and begR1 >= begR2 + distance) {
				candidates[reg2][reg1]++;
			}
		}
		i++;
	}
	countMtx.lock();
	processedBarcodes++;
	fprintf(stderr, " \rProcessed %d / %d barcodes (%d%%)", processedBarcodes, totalBarcodes, (int) ((float) processedBarcodes / totalBarcodes * 100));
	countMtx.unlock();

	return candidates;
}

pair<robin_hood::unordered_set<string>, vector<BamAlignment>> extractAlignmentsAndHeadersFromRegion(string& bamFile, string& region, unsigned additionalSize) { 
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.LocateIndex()) {
		fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	// Move to the region of interest, extending it if needed
	if (additionalSize > 0) {
		vector<string> t = splitString(region, ":");
	    vector<string> p = splitString(t[1], "-");
		string region1 = t[0] + ":" + to_string(max(static_cast<uint32_t>(0), static_cast<uint32_t>(stoul(p[0]) - additionalSize))) + "-" + to_string(static_cast<uint32_t>(stoul(p[1]) + additionalSize));
		reader.SetRegion(stringToBamRegion(reader, region1));	
	} else {
		reader.SetRegion(stringToBamRegion(reader, region));
	}
	

	robin_hood::unordered_set<string> headers;
	vector<BamAlignment> alignments;
	BamAlignment alignment;
	while (reader.GetNextAlignment(alignment)) {
		if (alignment.IsMapped()) {
			headers.insert(alignment.Name);
			alignments.push_back(alignment);
		}
	}

	return make_pair(headers, alignments);
}

vector<BamAlignment> retrieveAlignmentsWithCommonHeadersAndTreatInsertions(vector<BamAlignment>& alignments, robin_hood::unordered_set<string>& headers, SVSupports& support, unsigned posSupportRegion[], int32_t beg) {
	vector<BamAlignment> common;

	for (auto al : alignments) {
		if (headers.count(al.Name)) {
			common.push_back(al);
		} else {
			// TODO vérifier ça
			if (al.IsPrimaryAlignment() and !al.IsMateMapped()) {
				support.insertion++;
			} else if (!al.IsMateMapped()) {
				unsigned pos = al.Position;
				std::vector<CigarOp> cigar = al.CigarData;
				unsigned i = 0;
				while (i < cigar.size() and cigar[i].Type == 'M') {
					pos += cigar[i].Length;
					i++;
				}
				posSupportRegion[pos - beg] += 1;
				support.splitReads++;
			}
		}
	}

	return common;
}

bool orderBamAlignments (BamAlignment al1, BamAlignment al2) { 
	return al1.Name < al2.Name;
}

unsigned computeMaxSupport(SVSupports s) {
	return max(s.deletion, max(s.duplication, max(s.inversion, max(s.insertion, s.translocation))));
}

robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> processBarcodes(int nbThreads, robin_hood::unordered_map<string, int>& refIDs, vector<string>& regionsList, unsigned nbBins, string& bamFile, BarcodesPositionsIndex& barcodesPositionsIndex, int barcodesSize, int distance) {
	robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> candidates;
	unsigned poolSize = 10000;
	ctpl::thread_pool myPool(nbThreads);
	unsigned jobsLoaded = 0;
	unsigned jobsCompleted = 0;
	robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> tmpRes;
	robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> tmpCandidates;

	string begRegion = regionsList[0];
	vector<string> v1 = splitString(begRegion, ":");
	vector<string> vv1 = splitString(v1[1], "-");
	pair<int32_t, int32_t> beg = make_pair(refIDs[v1[0]], static_cast<uint32_t>(stoul(vv1[0])));
	string endRegion;
	for (unsigned i = 1; i <= nbBins; i++) {
		cerr << "Processing barcodes, iteration " << i << " of " << nbBins << endl;
		processedBarcodes = 0;
		endRegion = regionsList[i * regionsList.size() / nbBins - 1];
		vector<string> v2 = splitString(endRegion, ":");
		vector<string> vv2 = splitString(v2[1], "-");
		pair<int32_t, int32_t> end = make_pair(refIDs[v2[0]], static_cast<uint32_t>(stoul(vv2[0])));
		tmpCandidates.clear();

		jobsLoaded = 0;
		jobsCompleted = 0;

		auto it = barcodesPositionsIndex.begin();

		// Load the first jobs
		vector<std::future<robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>>> results(poolSize);
	    while (jobsLoaded < poolSize && jobsLoaded < barcodesPositionsIndex.size()) {
	        results[jobsLoaded] = myPool.push(processBarcode, ref(beg), ref(end), ref(bamFile), ref(barcodesPositionsIndex), ref(it->first), barcodesPositionsIndex.size(), ref(distance));
	        jobsLoaded++;
	        it++;
		}

		// Load the remaining jobs as other jobs terminate
		unsigned curJob = 0;
	    std::pair<std::string, std::string> curRes;
	    while(jobsLoaded < barcodesPositionsIndex.size()) {
	    	// Get the job results
	        tmpRes = results[curJob].get();
	        // Fill the results map
	        for (auto p : tmpRes) {
	        	for (auto pp : p.second) {
		        	tmpCandidates[p.first][pp.first] += pp.second;
		        }
	        }
	        jobsCompleted++;
	        
	        // Load the next job
	        results[curJob] = myPool.push(processBarcode, ref(beg), ref(end), ref(bamFile), ref(barcodesPositionsIndex), ref(it->first), barcodesPositionsIndex.size(), ref(distance));
	        jobsLoaded++;
	        it++;
	        
	        // Increment the current job nb, and loop if needed
	        curJob++;
	        if(curJob == poolSize) {
	            curJob = 0;
	        }
		}


		// Wait for the remaining jobs to terminate
		while(jobsCompleted < jobsLoaded) {
	        // Get the job results
	        tmpRes = results[curJob].get();
	        // Fill the results map
	        for (auto p : tmpRes) {
	        	for (auto pp : p.second) {
		        	tmpCandidates[p.first][pp.first] += pp.second;
		        }
	        }
	        jobsCompleted++;
	        
	        // Increment the current job nb, and loop if needed
	        curJob++;
	        if(curJob == poolSize) {
	            curJob = 0;
	        }
		}

		// Store this step's candidates SV regions if they have a support higer than one
		for (auto c : tmpCandidates) {
			for (auto cc : c.second) {
				if (cc.second > 1) {
					candidates[c.first][cc.first] += cc.second;
				}
			}
		}
		
		begRegion = endRegion;
		beg = end;

		cerr << endl;
	}

	return candidates;
}

SVSupports computePairedReadsSupport(vector<BamAlignment>& commonAlignments1, vector<BamAlignment>& commonAlignments2, unsigned begR1, unsigned begR2, unsigned posSupportRegion1[], unsigned posSupportRegion2[]) {
	SVSupports supports;
	unsigned i = 0;
	unsigned j = 0;
	unsigned oldJ = 0;

	supports.alignments1 = commonAlignments1.size();
	supports.alignments2 = commonAlignments2.size();

	vector<int> clipSizes;
	vector<int> readPositions;
	vector<int> genomePositions;

	// TODO ATTENTION ! On peut compter plusieurs fois les splits ici... Car on fait des paires

	while (i < commonAlignments1.size()) {
		oldJ = j;
		while (j < commonAlignments2.size() and commonAlignments1[i].Name == commonAlignments2[j].Name) {
			if (commonAlignments1[i].IsPrimaryAlignment() and commonAlignments2[j].IsPrimaryAlignment()) {
				if (commonAlignments1[i].RefID != commonAlignments2[j].RefID) { 
					supports.translocation++;
				} else {
					if (commonAlignments1[i].IsFirstMate() and !commonAlignments1[i].IsReverseStrand() and commonAlignments2[j].IsSecondMate() and commonAlignments2[j].IsReverseStrand() and commonAlignments1[i].Position < commonAlignments2[j].Position) {
						supports.deletion++;
					}
					if (commonAlignments2[j].IsFirstMate() and !commonAlignments2[j].IsReverseStrand() and commonAlignments1[i].IsSecondMate() and commonAlignments1[i].IsReverseStrand()  and commonAlignments2[j].Position < commonAlignments1[i].Position) {
						supports.deletion++;
					}
					if (commonAlignments1[i].IsFirstMate() and commonAlignments1[i].IsReverseStrand() and commonAlignments2[j].IsSecondMate() and !commonAlignments2[j].IsReverseStrand() and commonAlignments1[i].Position < commonAlignments2[j].Position) {
						supports.deletion++;
					}
					if (commonAlignments2[j].IsFirstMate() and commonAlignments2[j].IsReverseStrand() and commonAlignments1[i].IsSecondMate() and !commonAlignments1[i].IsReverseStrand()  and commonAlignments2[j].Position < commonAlignments1[i].Position) {
						supports.deletion++;
					}

					if (commonAlignments1[i].IsFirstMate() and !commonAlignments1[i].IsReverseStrand() and commonAlignments2[j].IsSecondMate() and commonAlignments2[j].IsReverseStrand() and commonAlignments2[j].Position < commonAlignments1[i].Position) {
						supports.duplication++;
					}
					if (commonAlignments2[j].IsFirstMate() and !commonAlignments2[j].IsReverseStrand() and commonAlignments1[i].IsSecondMate() and commonAlignments1[i].IsReverseStrand() and commonAlignments1[i].Position < commonAlignments2[j].Position) {
						supports.duplication++;
					}
					if (commonAlignments1[i].IsFirstMate() and commonAlignments1[i].IsReverseStrand() and commonAlignments2[j].IsSecondMate() and !commonAlignments2[j].IsReverseStrand() and commonAlignments2[j].Position < commonAlignments1[i].Position) {
						supports.duplication++;
					}
					if (commonAlignments2[j].IsFirstMate() and commonAlignments2[j].IsReverseStrand() and commonAlignments1[i].IsSecondMate() and !commonAlignments1[i].IsReverseStrand() and commonAlignments1[i].Position < commonAlignments2[j].Position) {
						supports.duplication++;
					}

					if (commonAlignments1[i].IsFirstMate() and !commonAlignments1[i].IsReverseStrand() and commonAlignments2[j].IsSecondMate() and !commonAlignments2[j].IsReverseStrand()) {
						supports.inversion++;
					}
					if (commonAlignments2[j].IsFirstMate() and !commonAlignments2[j].IsReverseStrand() and commonAlignments1[i].IsSecondMate() and !commonAlignments1[i].IsReverseStrand()) {
						supports.inversion++;
					}
					if (commonAlignments1[i].IsFirstMate() and commonAlignments1[i].IsReverseStrand() and commonAlignments2[j].IsSecondMate() and commonAlignments2[j].IsReverseStrand()) {
						supports.inversion++;
					}
					if (commonAlignments2[j].IsFirstMate() and commonAlignments2[j].IsReverseStrand() and commonAlignments1[i].IsSecondMate() and commonAlignments1[i].IsReverseStrand()) {
						supports.inversion++;
					}
				} 

			}

			// TODO prendre en compte les signature split pour le typage

			if (commonAlignments1[i].GetSoftClips(clipSizes, readPositions, genomePositions, false) and commonAlignments2[j].GetSoftClips(clipSizes, readPositions, genomePositions, false)) {
				// cerr << "YES !!!" << endl;
			}
			
			if (commonAlignments1[i].GetSoftClips(clipSizes, readPositions, genomePositions, false)) {
				BamAlignment a = commonAlignments1[i];
				unsigned pos = a.Position;
				std::vector<CigarOp> cigar = a.CigarData;
				unsigned i = 0;
				while (i < cigar.size() and cigar[i].Type == 'M') {
					pos += cigar[i].Length;
					i++;
				}
				posSupportRegion1[pos - begR1] += 1;
				supports.splitReads++;
			}			
			
			if (commonAlignments2[j].GetSoftClips(clipSizes, readPositions, genomePositions, false)) {
				BamAlignment a = commonAlignments2[j];
				unsigned pos = a.Position;
				std::vector<CigarOp>cigar = a.CigarData;
				unsigned i = 0;
				while (i < cigar.size() and cigar[i].Type == 'M') {
					pos += cigar[i].Length;
					i++;
				}
				posSupportRegion2[pos - begR2] += 1;
				supports.splitReads++;
			}

			j++;
		}

		if (i < commonAlignments1.size() - 1 and commonAlignments1[i].Name == commonAlignments1[i+1].Name) {
			j = oldJ;
		}
		i++;
	}

	return supports;
}

pair<robin_hood::unordered_set<string>, vector<BamAlignment>> extractAlignmentsWithCommonHeadersAndTreatInsertions(string& bamFile, string& region, unsigned additionalSize, robin_hood::unordered_set<string>& headers, SVSupports& support, unsigned posSupportRegion[], int32_t beg) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.LocateIndex()) {
		fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	// Move to the region of interest, extending it if needed
	if (additionalSize > 0) {
		vector<string> t = splitString(region, ":");
	    vector<string> p = splitString(t[1], "-");
		string region1 = t[0] + ":" + to_string(max(static_cast<uint32_t>(0), static_cast<uint32_t>(stoul(p[0]) - additionalSize))) + "-" + to_string(static_cast<uint32_t>(stoul(p[1]) + additionalSize));
		reader.SetRegion(stringToBamRegion(reader, region1));	
	} else {
		reader.SetRegion(stringToBamRegion(reader, region));
	}

	robin_hood::unordered_set<string> headersSet;
	vector<BamAlignment> alignments;
	BamAlignment alignment;
	while (reader.GetNextAlignment(alignment)) {
		if (alignment.IsMapped() and headers.count(alignment.Name) != 0) {
			headersSet.insert(alignment.Name);
			alignments.push_back(alignment);
		} else {
			// TODO vérifier ça (insertion)
			if (alignment.IsPrimaryAlignment() and !alignment.IsMateMapped()) {
				support.insertion++;
			} else if (!alignment.IsMateMapped()) {
				unsigned pos = alignment.Position;
				std::vector<CigarOp> cigar = alignment.CigarData;
				unsigned i = 0;
				while (i < cigar.size() and cigar[i].Type == 'M') {
					pos += cigar[i].Length;
					i++;
				}
				posSupportRegion[pos - beg] += 1;
				support.splitReads++;
			}
		}
	}

	return make_pair(headersSet, alignments);
}

vector<pair<pair<string*, string*>, SVSupports>> processCandidate(int id, pair<string*, robin_hood::unordered_map<string*, unsigned>>& candidate, string& bamFile) {
	vector<pair<pair<string*, string*>, SVSupports>> res;

	// Start positions and arrays for breakpoints computation
    int32_t begR1 = regionBegPosition(*(candidate.first)) - 2 * windowSize;

    unsigned posSupportRegion1[5 * windowSize];
	for (unsigned i = 0; i < 5 * windowSize; i++) {
		posSupportRegion1[i] = 0;
	}

	pair<robin_hood::unordered_set<string>, vector<BamAlignment>> hr = extractAlignmentsAndHeadersFromRegion(bamFile, *(candidate.first), 0);
	robin_hood::unordered_set<string> headers1 = hr.first;
	vector<BamAlignment> alignments1 = hr.second;

	for (auto p : candidate.second) {
		for (unsigned i = 0; i < 5 * windowSize; i++) {
			posSupportRegion1[i] = 0;
		}

		SVSupports support;

		int32_t begR2 = regionBegPosition(*(p.first)) - 2 * windowSize;

		unsigned posSupportRegion2[5 * windowSize];
		for (unsigned i = 0; i < 5 * windowSize; i++) {
			posSupportRegion2[i] = 0;
		}

		hr = extractAlignmentsWithCommonHeadersAndTreatInsertions(bamFile, *(p.first), 0, headers1, support, posSupportRegion2, begR2);
		robin_hood::unordered_set<string> headers2 = hr.first;
		vector<BamAlignment> commonAlignments2 = hr.second;

		vector<BamAlignment> commonAlignments1 = retrieveAlignmentsWithCommonHeadersAndTreatInsertions(alignments1, headers2, support, posSupportRegion1, begR1);

		std::sort(commonAlignments1.begin(), commonAlignments1.end(), orderBamAlignments);
		std::sort(commonAlignments2.begin(), commonAlignments2.end(), orderBamAlignments);
		
		support = computePairedReadsSupport(commonAlignments1, commonAlignments2, begR1, begR2, posSupportRegion1, posSupportRegion2);
		support.barcodes = p.second; 

		unsigned maxR1 = 0;
		unsigned j = 1;

		while (j < 5 * windowSize) {
			if (posSupportRegion1[j] > posSupportRegion1[maxR1]) {
				maxR1 = j;
			}
			j++;
		}

		unsigned maxR2 = 0;
		j = 1;

		while (j < 5 * windowSize) {
			if (posSupportRegion2[j] > posSupportRegion2[maxR2]) {
				maxR2 = j;
			}
			j++;
		}

		support.breakpoint1 = maxR1 + begR1;
		support.breakpoint2 = maxR2 + begR2;
		support.support1 = posSupportRegion1[maxR1];
		support.support2 = posSupportRegion2[maxR2];

		if (support.support1 > 0 and support.support2 > 0 and (support.translocation > 0 or support.deletion > 0 or support.duplication > 0 or support.inversion > 0 or support.insertion > 0)) {
			res.push_back(make_pair(make_pair(candidate.first, p.first), support));
		}

			countMtx.lock();
			processedCandidates++;
			fprintf(stderr, " \rProcessed %d / %d SV candidates (%d%%)", processedCandidates, totalCandidates, (int) ((float) processedCandidates / totalCandidates * 100));
			countMtx.unlock();
	}

	return res;
}

robin_hood::unordered_map<pair<string*, string*>, SVSupports, hashPointersPairs> processCandidates(unsigned nbThreads, robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates, string& bamFile) {
	robin_hood::unordered_map<pair<string*, string*>, SVSupports, hashPointersPairs> calledSVs;
	unsigned poolSize = 100000;
	ctpl::thread_pool myPool(nbThreads);
	unsigned jobsLoaded = 0;
	unsigned jobsCompleted = 0;
	vector<pair<pair<string*, string*>, SVSupports>> tmpRes;

	auto it = candidates.begin();

	// Load the first jobs
	vector<std::future<vector<pair<pair<string*, string*>, SVSupports>>>> results(poolSize);
    while (jobsLoaded < poolSize && jobsLoaded < candidates.size()) {
        results[jobsLoaded] = myPool.push(processCandidate, make_pair(it->first, it->second), ref(bamFile));
        jobsLoaded++;
        it++;
	}

	// Load the remaining jobs as other jobs terminate
	unsigned curJob = 0;
    std::pair<std::string, std::string> curRes;
    while(jobsLoaded < candidates.size()) {
    	// Get the job results
        tmpRes = results[curJob].get();
        // Fill the results map
        for (auto r : tmpRes) {
        	calledSVs[r.first] = r.second;
        }
        jobsCompleted++;
        
        // Load the next job
        results[curJob] = myPool.push(processCandidate, make_pair(it->first, it->second), ref(bamFile));
        jobsLoaded++;
        it++;
        
        // Increment the current job nb, and loop if needed
        curJob++;
        if(curJob == poolSize) {
            curJob = 0;
        }
	}


	// Wait for the remaining jobs to terminate
	while(jobsCompleted < jobsLoaded) {
        // Get the job results
        tmpRes = results[curJob].get();
        // Fill the results map
		for (auto r : tmpRes) {
        	calledSVs[r.first] = r.second;
        }
        jobsCompleted++;
        
        // Increment the current job nb, and loop if needed
        curJob++;
        if(curJob == poolSize) {
            curJob = 0;
        }
	}

	cerr << endl;

	return calledSVs;
}

bool compareSVSupport(const pair<pair<string*, string*>, SVSupports>& a, const pair<pair<string*, string*>, SVSupports>& b) { 
    // return a.second.barcodes + computeMaxSupport(a.second) > b.second.barcodes + computeMaxSupport(b.second); 
    return a.second.barcodes > b.second.barcodes; 
}

unsigned analyzeDistribution(robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates) {
	unsigned maxSupp = 0;
	for (auto p : candidates) {
		for (auto pp : p.second) {
			maxSupp = pp.second > maxSupp ? pp.second : maxSupp;
		}
	}
	maxSupp++;

	unsigned distrib[maxSupp];
	for (unsigned i = 0; i < maxSupp; i++) {
		distrib[i] = 0;
	}

	unsigned long totalSupp = 0;

	for (auto p : candidates) {
		for (auto pp : p.second) {
			for (unsigned j = pp.second; j <= pp.second; j++) {
				distrib[j]++;
				totalSupp++;
			}
		}
	}

	unsigned average = 1;
	while (average < maxSupp and distrib[average] == 0) {
		average++;
	}
	totalSupp = distrib[average];
	average++;
	while (average < maxSupp and (double) (totalSupp + distrib[average]) / totalSupp * 100 >= 100.001) {
		totalSupp += distrib[average];
		average++;
	}

	return average;
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
		if (reader.SetRegion(id, 0, id, d.RefLength - 1) and reader.GetNextAlignment(al)) {
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

void removeInvalidCandidates(robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates, unsigned average, unsigned maxRegionsLinks) {
	robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> filteredCandidates;

 	// Remove candidates that do not have enough support
	for (auto p : candidates) {
		for (auto pp : p.second) {
			if (pp.second >= average) {
				filteredCandidates[p.first][pp.first] += pp.second;
			}
		}
	}

	candidates.clear();

	// Remove regions that share similarities with too many other regions
	for (auto p : filteredCandidates) {
		if (p.second.size() < maxRegionsLinks) {
			candidates[p.first] = p.second;
			totalCandidates += p.second.size();
		}
	}
 }

robin_hood::unordered_set<string> validatesSVs(robin_hood::unordered_map<pair<string*, string*>, SVSupports, hashPointersPairs>& calledSVs) {
 	vector<pair<pair<string*, string*>, SVSupports>> sortedSVs;
	robin_hood::unordered_set<pair<string, unsigned>, hashPairs> consideredBreakpoints;
	robin_hood::unordered_set<string> finalSVs;

 	for (auto c : calledSVs) {
		sortedSVs.push_back(make_pair(c.first, c.second));
	}
	sort(sortedSVs.begin(), sortedSVs.end(), compareSVSupport);

	for (pair<pair<string*, string*>, SVSupports> sv : sortedSVs) {
		string s1 = splitString(*(sv.first.first), ":")[0];
		string s2 = splitString(*(sv.first.second), ":")[0];	
		if (consideredBreakpoints.count(make_pair(s1, sv.second.breakpoint1)) == 0 and consideredBreakpoints.count(make_pair(s2, sv.second.breakpoint2)) == 0) {
			string type = "DEL";
			unsigned maxSupp = sv.second.deletion;
			if (sv.second.duplication > maxSupp) {
				type = "DUP";
				maxSupp = sv.second.duplication;
			}
			if (sv.second.inversion > maxSupp) {
				type = "INV";
				maxSupp = sv.second.inversion;
			}
			if (sv.second.insertion > maxSupp) {
				type = "INS";
				maxSupp = sv.second.insertion;
			}
			if (sv.second.translocation > maxSupp) {
				type = "TRA";
				maxSupp = sv.second.translocation;
			}
			// double diffSupport = (double) max(sv.second.support1, sv.second.support2) / min(sv.second.support1, sv.second.support2);
			// if (maxSupp > 1 and sv.second.support1 > 1 and sv.second.support2 > 1) {
				string s = s1 + "\t" + to_string(sv.second.breakpoint1) + "\t" + s2 + "\t" + to_string(sv.second.breakpoint2) + "\t" + type;
				finalSVs.insert(s);
				
				// cerr << *(sv.first.first) << " " << *(sv.first.second) << endl;
				// cerr << "Alignments in region 1 : " << sv.second.alignments1 << endl;
				// cerr << "Alignments in region 2 : " << sv.second.alignments2 << endl;
				// cerr << "Barcodes : " << sv.second.barcodes << endl;
				// cerr << "Deletion : " << sv.second.deletion << endl;
				// cerr << "Duplication : " << sv.second.duplication << endl;
				// cerr << "Inversion : " << sv.second.inversion << endl;
				// cerr << "Insertion : " << sv.second.insertion << endl;
				// cerr << "Translocation : " << sv.second.translocation << endl;
				// cerr << "Breakpoint1 : " << sv.second.breakpoint1 << endl;
				// cerr << "Breakpoint2 : " << sv.second.breakpoint2 << endl;
				// cerr << "Support1 : " << sv.second.support1 << endl;
				// cerr << "Support2 : " << sv.second.support2 << endl;
				// cerr << sv.second.splitReads << endl;
				// cerr << endl;
				consideredBreakpoints.insert(make_pair(s1, sv.second.breakpoint1));
				consideredBreakpoints.insert(make_pair(s2, sv.second.breakpoint2));
			// }
		}
	}

	return finalSVs;
}

void outputSVs(robin_hood::unordered_set<string>& finalSVs) {
	for (string s : finalSVs) {
		cout << s << endl;
	}
}

int main(int argc, char* argv[]) {
	// Parse arguments
	string bamFile = argv[1];
	// string contig = argv[2];
	// string region = argv[3];
	string barcodesOffsetsIndexFile = argv[2];
	windowSize = atol(argv[3]);
	int distance = atoi(argv[4]);
	int nbThreads = atoi(argv[5]);
	vector<string> s1;
	vector<string> s2;
	unsigned nbBins = 10;
	unsigned maxRegionsLinks = 50;
	totalCandidates = 0;
	string validCandidatesFile = "candidates";

	cerr << "Preparing auxiliary data" << endl;
	prepareAuxiliaryData(bamFile);

	cerr << "Loading the barcodes positions index" << endl;
	BarcodesPositionsIndex barcodesPositionsIndex = loadBarcodesPositionsIndex(barcodesOffsetsIndexFile);
	// Remove adjacent positions that correspond to the same window
	removeIndexRedundancy(barcodesPositionsIndex);
	
	// Process every barcode to look for SV evidence
	totalBarcodes = barcodesPositionsIndex.size();
	robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> candidates = processBarcodes(nbThreads, refIDs, regionsList, nbBins, bamFile, barcodesPositionsIndex, barcodesPositionsIndex.size(), distance);

	cerr << "Computing and analyzing the distribution of shared barcodes between pairs of regions" << endl;
	unsigned average = analyzeDistribution(candidates);	

	cerr << "Removing invalid candidates (regions pairs that share less than " << average << " barcodes and regions that are paired with more than " << maxRegionsLinks << " other regions)" << endl;
	removeInvalidCandidates(candidates, average, maxRegionsLinks);

	cerr << "Saving SV candidates to file " << validCandidatesFile << endl;
	saveSVCandidates(validCandidatesFile, candidates);

	// Analyze alignments of valid candidates
	cerr << "Number of SV candidates to consider : " << totalCandidates << endl;
	robin_hood::unordered_map<pair<string*, string*>, SVSupports, hashPointersPairs> calledSVs = processCandidates(nbThreads, candidates, bamFile);
	candidates.clear();

	robin_hood::unordered_set<string> finalSVs = validatesSVs(calledSVs);

	cerr << "Output " << finalSVs.size() << " SVs" << endl;
	outputSVs(finalSVs);

	return EXIT_SUCCESS;
}

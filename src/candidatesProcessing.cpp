#include "candidatesProcessing.h"
#include "ctpl_stl.h"
#include <future>
#include "globalVariables.h"
#include "alignmentsProcessing.h"
#include "supportComputation.h"
#include <chrono>

using namespace std::chrono;

vector<pair<pair<string*, string*>, SVSupports>> processCandidate(int id, pair<string*, robin_hood::unordered_map<string*, unsigned>>& candidate, string& bamFile) {
	auto start = high_resolution_clock::now();

	vector<pair<pair<string*, string*>, SVSupports>> res;

	// Start positions and arrays for breakpoints computation
    int32_t begR1 = regionBegPosition(*(candidate.first)) - windowSize;

    unsigned posSupportRegion1[3 * windowSize];
	// for (unsigned i = 0; i < 3 * windowSize; i++) {
	// 	posSupportRegion1[i] = 0;
	// }

	pair<robin_hood::unordered_set<string>, vector<BamAlignment>> hr = extractAlignmentsAndHeadersFromRegion(bamFile, *(candidate.first), 0);
	robin_hood::unordered_set<string> headers1 = hr.first;
	vector<BamAlignment> alignments1 = hr.second;

	for (auto p : candidate.second) {
		for (unsigned i = 0; i < 3 * windowSize; i++) {
			posSupportRegion1[i] = 0;
		}

		SVSupports support;

		int32_t begR2 = regionBegPosition(*(p.first)) - windowSize;

		unsigned posSupportRegion2[3 * windowSize];
		for (unsigned i = 0; i < 3 * windowSize; i++) {
			posSupportRegion2[i] = 0;
		}

		hr = extractAlignmentsWithCommonHeadersAndTreatInsertions(bamFile, *(p.first), 0, headers1, support, posSupportRegion2, begR2);
		robin_hood::unordered_set<string> headers2 = hr.first;
		vector<BamAlignment> commonAlignments2 = hr.second;

		vector<BamAlignment> commonAlignments1 = retrieveAlignmentsWithCommonHeadersAndTreatInsertions(alignments1, headers2, support, posSupportRegion1, begR1);

		std::sort(commonAlignments1.begin(), commonAlignments1.end(), orderBamAlignments);
		std::sort(commonAlignments2.begin(), commonAlignments2.end(), orderBamAlignments);
		
		support = computeSVSupport(commonAlignments1, commonAlignments2, begR1, begR2, posSupportRegion1, posSupportRegion2);
		support.barcodes = p.second; 

		unsigned maxR1 = 0;
		unsigned j = 1;

		while (j < 3 * windowSize) {
			if (posSupportRegion1[j] > posSupportRegion1[maxR1]) {
				maxR1 = j;
			}
			j++;
		}

		unsigned maxR2 = 0;
		j = 1;

		while (j < 3 * windowSize) {
			if (posSupportRegion2[j] > posSupportRegion2[maxR2]) {
				maxR2 = j;
			}
			j++;
		}

		support.breakpoint1 = maxR1 + begR1;
		support.breakpoint2 = maxR2 + begR2;
		support.support1 = posSupportRegion1[maxR1];
		support.support2 = posSupportRegion2[maxR2];

		// commonAlignments1 = extractAlignmentsAndHeadersFromRegion(bamFile, *(candidate.first), 250).second;
		// commonAlignments2 = extractAlignmentsAndHeadersFromRegion(bamFile, *(p.first), 250).second;
		
		// // TODO : l'assemblage ne donne rien au dessus de k = 15
		// // if (support.support1 > 0 and support.support2 > 0) {
		// if (support.support1 > 0 and support.support2 > 0 and (support.translocation > 0 or support.deletion > 0 or support.duplication > 0 or support.inversion > 0 or support.insertion > 0)) {
		// 	string assembly; 
		// 	if (begR1 < begR2) {
		// 		assembly = computeAssembly(commonAlignments1, commonAlignments2, 31, 10, 100000, 5, *(candidate.first), *(p.first));
		// 	} else {
		// 		assembly = computeAssembly(commonAlignments2, commonAlignments1, 31, 10, 100000, 5, *(p.first), *(candidate.first));
		// 	}

		// 	cerr << "assembled : " << assembly << endl;
		// 	cerr << assembly.size() << endl;

		// 	pair<unsigned, unsigned> assemblyBreakpoints;
		// 	if (!assembly.empty()) {
		// 		assemblyBreakpoints = computeAssemblyBreakpoints(assembly, *(candidate.first), *(p.first));
		// 	}


		// 	cerr << "bp1 : " << support.breakpoint1 << " ; " << assemblyBreakpoints.first << endl;
		// 	cerr << "bp2 : " << support.breakpoint2 << " ; " << assemblyBreakpoints.second << endl;
		// }

		if (support.support1 > 0 and support.support2 > 0 and (support.translocation > 0 or support.deletion > 0 or support.duplication > 0 or support.inversion > 0 or support.insertion > 0)) {
			res.push_back(make_pair(make_pair(candidate.first, p.first), support));
		}

			countMtx.lock();
			processedCandidates++;
			fprintf(stderr, " \rProcessed %d / %d SV candidates (%d%%)", processedCandidates, totalCandidates, (int) ((float) processedCandidates / totalCandidates * 100));
			countMtx.unlock();
	}

	auto end = high_resolution_clock::now();

	auto duration = duration_cast<milliseconds>(end - start);
	// cerr << "Processing : " << *(candidate.first) << " (" << candidate.second.size() << " links) took : " << duration.count() << " ms" << endl;

	return res;
}

robin_hood::unordered_map<pair<string*, string*>, SVSupports, hashPointersPairs> processCandidates(unsigned nbThreads, robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates, string& bamFile) {
	robin_hood::unordered_map<pair<string*, string*>, SVSupports, hashPointersPairs> calledSVs;
	// nbThreads = 1;
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

void removeInvalidCandidates(robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates, Thresholds th, unsigned maxRegionsLinks) {
	robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> filteredCandidates;
	robin_hood::unordered_map<string*, unsigned> regionsLinks;

	unsigned small = 0;
	unsigned medium = 0;
	unsigned large = 0;

 	// Remove candidates that do not have enough support
	for (auto p : candidates) {
    	pair<string, int32_t> chrBeg1 = regionChrAndBegPosition(*(p.first));
    	string chr1 = chrBeg1.first;
    	int32_t beg1 = chrBeg1.second;
		for (auto pp : p.second) {
    		pair<string, int32_t> chrBeg2 = regionChrAndBegPosition(*(pp.first));
	    	string chr2 = chrBeg2.first;
    		int32_t beg2 = chrBeg2.second;

    		if (chr1 != chr2 and pp.second >= th.farTh) {
				filteredCandidates[p.first][pp.first] += pp.second;
				regionsLinks[p.first]++;
				regionsLinks[pp.first]++;
				large++;
        	} else if (chr1 == chr2 and abs(beg2 - beg1) <= mediumVariantsSize and pp.second >= th.closeTh) {
				filteredCandidates[p.first][pp.first] += pp.second;
				regionsLinks[p.first]++;
				regionsLinks[pp.first]++;
				small++;
        	} else if (chr1 == chr2 and abs(beg2 - beg1) > mediumVariantsSize and abs(beg2 - beg1) <= largeVariantsSize and pp.second >= th.averageTh) {
        		filteredCandidates[p.first][pp.first] += pp.second;
        		regionsLinks[p.first]++;
				regionsLinks[pp.first]++;
				medium++;
        	} else if (chr1 == chr2 and abs(beg2 - beg1) > largeVariantsSize and pp.second >= th.farTh) {
        		filteredCandidates[p.first][pp.first] += pp.second;
        		regionsLinks[p.first]++;
				regionsLinks[pp.first]++;
				large++;
        	}
		}
	}

	cerr << "small : " << small << endl;
	cerr << "medium : " << medium << endl;
	cerr << "large : " << large << endl;

	// candidates = filteredCandidates;
	candidates.clear();

	// Remove regions that share similarities with too many other regions
	unsigned addedLinks;
	for (auto p : filteredCandidates) {
		addedLinks = 0;
		if (regionsLinks[p.first] < maxRegionsLinks) {
			robin_hood::unordered_map<string*, unsigned> l;
			for (auto pp : p.second) {
				if (regionsLinks[pp.first] < maxRegionsLinks) {
					if (regionsLinks[p.first] > regionsLinks[pp.first]) {
						candidates[p.first][pp.first] = pp.second;
					} else {
						candidates[pp.first][p.first] = pp.second;
					}
					totalCandidates++;
				}
			}
		}
	}

 }

 void removeCandidatesRegionsLinks(robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates, unsigned maxRegionsLinks) {
	robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> filteredCandidates;
	robin_hood::unordered_map<string*, unsigned> regionsLinks;

	unsigned small = 0;
	unsigned medium = 0;
	unsigned large = 0;

	// Remove regions that share similarities with too many other regions
	unsigned addedLinks;
	for (auto p : candidates) {
		addedLinks = 0;
		if (regionsLinks[p.first] < maxRegionsLinks) {
			robin_hood::unordered_map<string*, unsigned> l;
			for (auto pp : p.second) {
				if (regionsLinks[pp.first] < maxRegionsLinks) {
					if (regionsLinks[p.first] > regionsLinks[pp.first]) {
						filteredCandidates[p.first][pp.first] = pp.second;
					} else {
						filteredCandidates[pp.first][p.first] = pp.second;
					}
					totalCandidates++;
				}
			}
		}
	}

	candidates = filteredCandidates;
 }

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

robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> loadSVCandidates(string file) {
	ifstream in;
	in.open(file);
	if (!in.is_open()) {
		fprintf(stderr, "Unable to open barcodes index file %s. Please provide an existing and valid file.\n", file.c_str());
		exit(EXIT_FAILURE);
	}

	robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> res;
	vector<string> v;
	string line;
	string* r1;
	string* r2;
	unsigned nbBarcodes;

	while (getline(in, line)) {
		v = splitString(line, "\t");
		// cerr << v[0] << " ; " << v[1] << " ; " << v[2] << endl;
		r1 = positionToRegion(windows[refIDs[v[0]]], stoul(v[1]));
		// cerr << *r1 << endl << endl;
		r2 = positionToRegion(windows[refIDs[v[3]]], stoul(v[4]));
		nbBarcodes = stoul(v[6]);
		res[r1][r2] = nbBarcodes;
		totalCandidates++;
	}

	in.close();
	return res;
}

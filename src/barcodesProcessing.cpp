#include "barcodesProcessing.h"
#include "ctpl_stl.h"
#include <future>
#include "globalVariables.h"
#include "misc.h"

unsigned long pairsWithOneBarcode;
unsigned long pairsWithMoreBarcodes;

// TODO TOUT PASSER EN REF AUX THREADS
robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> processBarcode(int id, pair<int32_t, int32_t> beg, pair<int32_t, int32_t> end, string& bamFile, BarcodesPositionsIndex& barcodesPositionsIndex, barcode& b, int barcodesSize, int& minVariantSize) {
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
			// TODO: Vérifier la minVariantSize entre les fenêtres
			// minVariantSize was 2 * windowSize
				// TODO inutile car beg1 sera toujours plus petit vu que l'index est trié ?
			} else if (reg1 != reg2 and begR2 >= begR1 + minVariantSize) {
				candidates[reg1][reg2]++;
			} else if (reg1 != reg2 and begR1 >= begR2 + minVariantSize) {
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

robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> processBarcodes(int nbThreads, robin_hood::unordered_map<string, int>& refIDs, vector<string>& regionsList, unsigned nbBins, string& bamFile, BarcodesPositionsIndex& barcodesPositionsIndex, int barcodesSize, int minVariantSize) {
	robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> candidates;
	ctpl::thread_pool myPool(nbThreads);
	unsigned jobsLoaded = 0;
	unsigned jobsCompleted = 0;
	robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> tmpRes;
	robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> tmpCandidates;
	pairsWithOneBarcode = 0;

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
	        results[jobsLoaded] = myPool.push(processBarcode, ref(beg), ref(end), ref(bamFile), ref(barcodesPositionsIndex), ref(it->first), barcodesPositionsIndex.size(), ref(minVariantSize));
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
	        results[curJob] = myPool.push(processBarcode, ref(beg), ref(end), ref(bamFile), ref(barcodesPositionsIndex), ref(it->first), barcodesPositionsIndex.size(), ref(minVariantSize));
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
				if (cc.second > minBarcodes) {
					candidates[c.first][cc.first] += cc.second;
					pairsWithMoreBarcodes++;
				} else {
					pairsWithOneBarcode++;
				}
			}
		}
		
		begRegion = endRegion;
		beg = end;

		cerr << endl;
	}

	cerr << "pairs with one barcode : " << pairsWithOneBarcode << endl;
	cerr << "pairs with more barcodes : " << pairsWithMoreBarcodes << endl;

	return candidates;
}

Thresholds analyzeDistribution(robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>>& candidates) {
    unsigned maxSupp = 0;
    for (auto p : candidates) {
        for (auto pp : p.second) {
            maxSupp = pp.second > maxSupp ? pp.second : maxSupp;
        }
    }
    maxSupp++;

    unsigned* closeDist = new unsigned[maxSupp];
    for (unsigned i = 0; i < maxSupp; i++) {
        closeDist[i] = 0;
    }

    unsigned* averageDist = new unsigned[maxSupp];
    for (unsigned i = 0; i < maxSupp; i++) {
        averageDist[i] = 0;
    }

    unsigned* farDist = new unsigned[maxSupp];
    for (unsigned i = 0; i < maxSupp; i++) {
        farDist[i] = 0;
    }

    unsigned long totalSuppClose = 0;
    unsigned long totalSuppAverage = 0;
    unsigned long totalSuppFar = 0;

    // TODO : Bug avec les seuils 99.95 ?

    for (auto p : candidates) {
    	pair<string, int32_t> chrBeg1 = regionChrAndBegPosition(*(p.first));
    	string chr1 = chrBeg1.first;
    	int32_t beg1 = chrBeg1.second;
        for (auto pp : p.second) {
	    	pair<string, int32_t> chrBeg2 = regionChrAndBegPosition(*(pp.first));
	    	string chr2 = chrBeg2.first;
    		int32_t beg2 = chrBeg2.second;
            for (unsigned j = pp.second; j <= pp.second; j++) {
            	if (chr1 != chr2) {
            		farDist[j]++;
            		totalSuppFar++;
            	} else if (abs(beg2 - beg1) <= mediumVariantsSize) {
            		closeDist[j]++;
            		totalSuppClose++;
            	} else if (abs(beg2 - beg1) <= largeVariantsSize) {
            		averageDist[j]++;
            		totalSuppAverage++;
            	} else {
            		farDist[j]++;
            		totalSuppFar++;
            	}
            }
        }
    }

    unsigned closeTh = 1;
    while (closeTh < maxSupp and closeDist[closeTh] == 0) {
        closeTh++;
    }
    unsigned long totalSupp = closeDist[closeTh];
    closeTh++;
    while (closeTh < maxSupp and ((double) totalSupp / totalSuppClose) * 100 < smallVariantsRate) {
    	totalSupp += closeDist[closeTh];
    	closeTh++;
    }

    unsigned averageTh = 1;
    while (averageTh < maxSupp and averageDist[averageTh] == 0) {
        averageTh++;
    }
    totalSupp = averageDist[averageTh];
    averageTh++;
    while (averageTh < maxSupp and ((double) totalSupp / totalSuppAverage) * 100 < mediumVariantsRate) {
    	totalSupp += averageDist[averageTh];
    	averageTh++;
    }

    unsigned farTh = 1;
    while (farTh < maxSupp and farDist[farTh] == 0) {
        farTh++;
    }
    totalSupp = farDist[farTh];
    farTh++;
    while (farTh < maxSupp and ((double) totalSupp / totalSuppFar) * 100 < largeVariantsRate) {
    	totalSupp += farDist[farTh];
    	farTh++;
    }

    Thresholds th;
    th.closeTh = closeTh;
    th.averageTh = averageTh;
    th.farTh = farTh;

    delete[] closeDist;
    delete[] averageDist;
    delete[] farDist;

    return th;
}
#include "SVValidation.h"
#include "misc.h"
#include "globalVariables.h"

bool isPresent(robin_hood::unordered_set<StructuralVariant>& SVs, StructuralVariant sv, unsigned distance) {
	for (StructuralVariant s : SVs) {


		if (s.chr1 == sv.chr1 and abs((int) s.breakpoint1 - (int) sv.breakpoint1) <= distance and s.chr2 == sv.chr2 and abs((int) s.breakpoint2 - (int) sv.breakpoint2) <= distance) {
			return true;
		}
	}

	return false;
}

bool breakpointIsPresent(robin_hood::unordered_set<StructuralVariant>& SVs, string chr, unsigned breakpoint, unsigned distance) {
	for (StructuralVariant s : SVs) {
		if ((s.chr1 == chr and s.breakpoint1 == breakpoint) or (s.chr2 == chr and s.breakpoint2 == breakpoint)) {
			return true;
		}
	}

	return false;
}

robin_hood::unordered_set<StructuralVariant> validatesSVs(robin_hood::unordered_map<pair<string*, string*>, SVSupports, hashPointersPairs>& calledSVs) {
 	vector<pair<pair<string*, string*>, SVSupports>> sortedSVs;
	robin_hood::unordered_set<pair<string, unsigned>, hashPairs> consideredBreakpoints;
	robin_hood::unordered_set<StructuralVariant> finalSVs;

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
				StructuralVariant s;
				if (s1 < s2) {
					s = StructuralVariant(s1, sv.second.breakpoint1, s2, sv.second.breakpoint2, type, sv.second.barcodes, maxSupp);
					// s = s1 + "\t" + to_string(sv.second.breakpoint1) + "\t" + s2 + "\t" + to_string(sv.second.breakpoint2) + "\t" + type;
				} else if (s2 < s1) {
					s = StructuralVariant(s2, sv.second.breakpoint2, s1, sv.second.breakpoint1, type, sv.second.barcodes, maxSupp);
					// s = s2 + "\t" + to_string(sv.second.breakpoint2) + "\t" + s1 + "\t" + to_string(sv.second.breakpoint1) + "\t" + type;
				} else if (sv.second.breakpoint1 < sv.second.breakpoint2) {
					s = StructuralVariant(s1, sv.second.breakpoint1, s2, sv.second.breakpoint2, type, sv.second.barcodes, maxSupp);
					// s = s1 + "\t" + to_string(sv.second.breakpoint1) + "\t" + s2 + "\t" + to_string(sv.second.breakpoint2) + "\t" + type;
				} else {
					s = StructuralVariant(s2, sv.second.breakpoint2, s1, sv.second.breakpoint1, type, sv.second.barcodes, maxSupp);
					// s = s2 + "\t" + to_string(sv.second.breakpoint2) + "\t" + s1 + "\t" + to_string(sv.second.breakpoint1) + "\t" + type;
				}
				if (!isPresent(finalSVs, s, duplicatesDistance)) {
					finalSVs.insert(s);
				
					cerr << *(sv.first.first) << " " << *(sv.first.second) << endl;
					cerr << "Alignments in region 1 : " << sv.second.alignments1 << endl;
					cerr << "Alignments in region 2 : " << sv.second.alignments2 << endl;
					cerr << "Barcodes : " << sv.second.barcodes << endl;
					cerr << "Deletion : " << sv.second.deletion << endl;
					cerr << "Duplication : " << sv.second.duplication << endl;
					cerr << "Inversion : " << sv.second.inversion << endl;
					cerr << "Insertion : " << sv.second.insertion << endl;
					cerr << "Translocation : " << sv.second.translocation << endl;
					cerr << "Breakpoint1 : " << sv.second.breakpoint1 << endl;
					cerr << "Breakpoint2 : " << sv.second.breakpoint2 << endl;
					cerr << "Support1 : " << sv.second.support1 << endl;
					cerr << "Support2 : " << sv.second.support2 << endl;
					cerr << sv.second.splitReads << endl;
					cerr << endl;
					consideredBreakpoints.insert(make_pair(s1, sv.second.breakpoint1));
					consideredBreakpoints.insert(make_pair(s2, sv.second.breakpoint2));
				}
			// }
		}
	}

	return finalSVs;
}


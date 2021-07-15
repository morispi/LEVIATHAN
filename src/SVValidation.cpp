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
			// if (maxSupp > 1 and sv.second.support1 > 1 and sv.second.support2 > 1) {
				StructuralVariant s;
				if (s1 < s2) {
					s = StructuralVariant(s1, sv.second.breakpoint1, s2, sv.second.breakpoint2, type, sv.second.barcodes, maxSupp);
				} else if (s2 < s1) {
					s = StructuralVariant(s2, sv.second.breakpoint2, s1, sv.second.breakpoint1, type, sv.second.barcodes, maxSupp);
				} else if (sv.second.breakpoint1 < sv.second.breakpoint2) {
					s = StructuralVariant(s1, sv.second.breakpoint1, s2, sv.second.breakpoint2, type, sv.second.barcodes, maxSupp);
				} else {
					s = StructuralVariant(s2, sv.second.breakpoint2, s1, sv.second.breakpoint1, type, sv.second.barcodes, maxSupp);
				}
				if (!isPresent(finalSVs, s, duplicatesDistance)) {
					finalSVs.insert(s);
					consideredBreakpoints.insert(make_pair(s1, sv.second.breakpoint1));
					consideredBreakpoints.insert(make_pair(s2, sv.second.breakpoint2));
				}
			// }
		}
	}

	return finalSVs;
}


#include "alignmentsProcessing.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "misc.h"

using namespace BamTools;

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
	vector<int> clipSizes;
	vector<int> readPositions;
	vector<int> genomePositions;

	for (auto al : alignments) {
		if (headers.count(al.Name)) {
			common.push_back(al);
		} else {
			// TODO vérifier ça
			// if (al.IsPrimaryAlignment() and !al.IsMateMapped()) {
			// 	support.insertion++;
			// } else if (!al.IsMateMapped()) {
			// 	unsigned pos = al.Position;
			// 	std::vector<CigarOp> cigar = al.CigarData;
			// 	unsigned i = 0;
			// 	while (i < cigar.size() and cigar[i].Type == 'M') {
			// 		pos += cigar[i].Length;
			// 		i++;
			// 	}
			// 	posSupportRegion[pos - beg] += 1;
			// 	support.splitReads++;
			// }
			if (!al.IsMateMapped()) {
				support.insertion++;
				if (al.GetSoftClips(clipSizes, readPositions, genomePositions, false)) {
					unsigned pos = al.Position;
					std::vector<CigarOp> cigar = al.CigarData;
					unsigned i = 0;
					while (i < cigar.size() and cigar[i].Type != 'S') {
						pos += cigar[i].Length;
						i++;
					}
					posSupportRegion[pos - beg] += 1;
					support.splitReads++;
				}
			}
		}
	}

	return common;
}

pair<robin_hood::unordered_set<string>, vector<BamAlignment>> extractAlignmentsWithCommonHeadersAndTreatInsertions(string& bamFile, string& region, unsigned additionalSize, robin_hood::unordered_set<string>& headers, SVSupports& support, unsigned posSupportRegion[], int32_t beg) {
	vector<int> clipSizes;
	vector<int> readPositions;
	vector<int> genomePositions;

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
	BamAlignment al;
	while (reader.GetNextAlignment(al)) {
		if (al.IsMapped() and headers.count(al.Name) != 0) {
			headersSet.insert(al.Name);
			alignments.push_back(al);
		} else {
			// TODO vérifier ça (insertion)
			// if (alignment.IsPrimaryAlignment() and !alignment.IsMateMapped()) {
			// 	support.insertion++;
			// } else if (!alignment.IsMateMapped()) {
			// 	unsigned pos = alignment.Position;
			// 	std::vector<CigarOp> cigar = alignment.CigarData;
			// 	unsigned i = 0;
			// 	while (i < cigar.size() and cigar[i].Type == 'M') {
			// 		pos += cigar[i].Length;
			// 		i++;
			// 	}
			// 	posSupportRegion[pos - beg] += 1;
			// 	support.splitReads++;
			// }
			if (!al.IsMateMapped()) {
				support.insertion++;
				if (al.GetSoftClips(clipSizes, readPositions, genomePositions, false)) {
					unsigned pos = al.Position;
					std::vector<CigarOp> cigar = al.CigarData;
					unsigned i = 0;
					while (i < cigar.size() and cigar[i].Type != 'S') {
						pos += cigar[i].Length;
						i++;
					}
					posSupportRegion[pos - beg] += 1;
					support.splitReads++;
				}
			}
		}
	}

	return make_pair(headersSet, alignments);
}
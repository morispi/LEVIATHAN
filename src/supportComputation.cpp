#include "supportComputation.h"
#include "api/BamAlignment.h"

using namespace BamTools;

SVSupports computeSVSupport(vector<BamAlignment>& commonAlignments1, vector<BamAlignment>& commonAlignments2, unsigned begR1, unsigned begR2, unsigned posSupportRegion1[], unsigned posSupportRegion2[]) {
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
				while (i < cigar.size() and cigar[i].Type != 'S') {
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
				while (i < cigar.size() and cigar[i].Type != 'S') {
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
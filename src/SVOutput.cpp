#include "SVOutput.h"
#include "robin_hood.h"
#include "genomeProcessing.h"
#include "globalVariables.h"

void outputSVsAsBed(robin_hood::unordered_set<StructuralVariant>& finalSVs) {
	for (StructuralVariant s : finalSVs) {
		cout << s.chr1 << "\t" << s.breakpoint1 << "\t" << s.chr2 << "\t" << s.breakpoint2 << "\t" << s.type << "\t" << endl;
	}
}

void outputSVsAsVCF(robin_hood::unordered_set<StructuralVariant>& finalSVs, string cmdLine, string reference, string outputFile) {
	ofstream out;
	out.open(outputFile, ios::out | ios::binary);
	if (!out.is_open()) {
		fprintf(stderr, "Unable to open file %s.", outputFile.c_str());
		exit(EXIT_FAILURE);
	}

	out << "##fileformat=VCFv4.2" << endl;
	out << "##source=MyVC" << endl;
	out << "##command=" << cmdLine << endl;
	out << "#reference=" << reference << endl;
	out << "##FILTER=<ID=PASS,Description=\"Filters passed\">" << endl;
	out << "##FILTER=<ID=LowQual,Description=\"Low quality\">" << endl;
	out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of variant, either DEL, DUP, INV, INS, or BND\">" << endl;
	out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl;
	out << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
	out << "##INFO=<ID=BARCODES,Number=1,Type=Integer,Description=\"The number of barcodes shared by the regions spanning the breakpoints of the variant\">" << endl;
	out << "##INFO=<ID=PAIRS,Number=1,Type=Integer,Description=\"The number of discordant read pairs validating the variant and its type\">" << endl;
	out << "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">" << endl;
	out << "#CHROM" << "\t" << "POS" << "\t" << "ID" << "\t" << "REF" << "\t" << "ALT" << "\t" << "QUAL" << "\t" << "FILTER" << "\t" << "INFO" << endl;
	
	unsigned calls = 0;
	for (StructuralVariant s : finalSVs) {
		if (s.type != "TRA") {
			out << s.chr1 << "\t" << s.breakpoint1 << "\t" << "call_" << calls << "\t" << twoBitsToString(genomeIndex[s.chr1][s.breakpoint1 * 2], genomeIndex[s.chr1][s.breakpoint1 * 2 + 1])
				 << "\t" << "<" << s.type << ">" << "\t" << s.barcodes + s.pairSupport << "\t" << "PASS" << "\t" << "SVTYPE=" << s.type << ";END=" << s.breakpoint2 
				 << ";SVLEN=" << s.breakpoint2 - s.breakpoint1 + 1 << ";BARCODES=" << s.barcodes << ";PAIRS=" << s.pairSupport << endl;
		} else {
			// TODO: gather 4 breakpoints for translocations?
			out << s.chr1 << "\t" << s.breakpoint1 << "\t" << "call_" << calls << "_1" << "\t" << twoBitsToString(genomeIndex[s.chr1][s.breakpoint1 * 2], genomeIndex[s.chr1][s.breakpoint1 * 2 + 1])
				 << "\t" << "]" << s.chr2 << ":" << s.breakpoint2 << "]" << twoBitsToString(genomeIndex[s.chr2][s.breakpoint2 * 2], genomeIndex[s.chr2][s.breakpoint2 * 2 + 1]) << "\t" 
				 << s.barcodes + s.pairSupport << "\t" << "PASS" << "\t" << "SVTYPE=BND" << ";MATEID=" << "call_" << calls << "_2" << ";BARCODES=" << s.barcodes << ";PAIRS=" << s.pairSupport << endl;

			out << s.chr2 << "\t" << s.breakpoint2 << "\t" << "call_" << calls << "_2" << "\t" << twoBitsToString(genomeIndex[s.chr2][s.breakpoint2 * 2], genomeIndex[s.chr2][s.breakpoint2 * 2 + 1])
				 << "\t" << twoBitsToString(genomeIndex[s.chr1][s.breakpoint1 * 2], genomeIndex[s.chr1][s.breakpoint1 * 2 + 1]) << "[" << s.chr1 << ":" << s.breakpoint1 << "[" 
				 << "\t" << s.barcodes + s.pairSupport << "\t" << "PASS" << "\t" << "SVTYPE=BND" << ";MATEID=" << "call_" << calls << "_1" << ";BARCODES=" << s.barcodes << ";PAIRS=" << s.pairSupport << endl;
		}

		calls++;
	}

}

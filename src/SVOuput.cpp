#include "SVOutput.h"
#include "robin_hood.h"
#include "genomeProcessing.h"
#include "globalVariables.h"

void outputSVsAsBed(robin_hood::unordered_set<StructuralVariant>& finalSVs) {
	for (StructuralVariant s : finalSVs) {
		cout << s.chr1 << "\t" << s.breakpoint1 << "\t" << s.chr2 << "\t" << s.breakpoint2 << "\t" << s.type << "\t" << endl;
	}
}

void outputSVsAsVCF(robin_hood::unordered_set<StructuralVariant>& finalSVs, string cmdLine, string reference) {
	cout << "##fileformat=VCFv4.2" << endl;
	cout << "##source=MyVC" << endl;
	cout << "##command=" << cmdLine << endl;
	cout << "#reference=" << reference << endl;
	cout << "##FILTER=<ID=PASS,Description=\"Filters passed\">" << endl;
	cout << "##FILTER=<ID=LowQual,Description=\"Low quality\">" << endl;
	cout << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of variant, either DEL, DUP, INV, INS, or BND\">" << endl;
	cout << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl;
	cout << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
	cout << "##INFO=<ID=BARCODES,Number=1,Type=Integer,Description=\"The number of barcodes shared by the regions spanning the breakpoints of the variant\">" << endl;
	cout << "##INFO=<ID=PAIRS,Number=1,Type=Integer,Description=\"The number of discordant read pairs validating the variant and its type\">" << endl;
	cout << "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">" << endl;
	cout << "#CHROM" << "\t" << "POS" << "\t" << "ID" << "\t" << "REF" << "\t" << "ALT" << "\t" << "QUAL" << "\t" << "FILTER" << "\t" << "INFO" << endl;
	
	unsigned calls = 0;
	for (StructuralVariant s : finalSVs) {
		if (s.type != "TRA") {
			cout << s.chr1 << "\t" << s.breakpoint1 << "\t" << "call_" << calls << "\t" << fullnum2str(genomeIndex[s.chr1]).substr(s.breakpoint1, 1) 
				 << "\t" << "<" << s.type << ">" << "\t" << s.barcodes + s.pairSupport << "\t" << "PASS" << "\t" << "SVTYPE=" << s.type << ";END=" << s.breakpoint2 
				 << ";SVLEN=" << s.breakpoint2 - s.breakpoint1 + 1 << ";" << ";BARCODES=" << s.barcodes << ";PAIRS=" << s.pairSupport << endl;
		} else {
			cout << s.chr1 << "\t" << s.breakpoint1 << "\t" << "call_" << calls << "_1" << "\t" << fullnum2str(genomeIndex[s.chr1]).substr(s.breakpoint1, 1) 
				 << "\t" << "]" << s.chr2 << ":" << s.breakpoint2 << "]" << fullnum2str(genomeIndex[s.chr2]).substr(s.breakpoint2, 1) << ">" << "\t" 
				 << s.barcodes + s.pairSupport << "\t" << "PASS" << "\t" << "SVTYPE=BND" << ";MATEID=" << "call_" << calls << "_2" << ";END=" << s.breakpoint2 
				 << ";SVLEN=" << s.breakpoint2 - s.breakpoint1 + 1 << ";" << ";BARCODES=" << s.barcodes << ";PAIRS=" << s.pairSupport << endl;

			cout << s.chr1 << "\t" << s.breakpoint1 << "\t" << "call_" << calls << "_2" << "\t" << fullnum2str(genomeIndex[s.chr1]).substr(s.breakpoint1, 1) 
				 << "\t" << fullnum2str(genomeIndex[s.chr2]).substr(s.breakpoint2, 1) << "[" << s.chr2 << ":" << s.breakpoint2 << "[" 
				 << s.barcodes + s.pairSupport << "\t" << "PASS" << "\t" << "SVTYPE=BND" << ";MATEID=" << "call_" << calls << "_1" << ";END=" << s.breakpoint2 
				 << ";SVLEN=" << s.breakpoint2 - s.breakpoint1 + 1 << ";" << ";BARCODES=" << s.barcodes << ";PAIRS=" << s.pairSupport << endl;
		}

		calls++;
	}

}
#include "supportComputation.h"
#include "alignmentsProcessing.h"
#include "candidatesProcessing.h"
#include "barcodesProcessing.h"
#include "misc.h"
#include "globalVariables.h"
#include "SVSupports.h"
#include "genomeProcessing.h"
#include "SVOutput.h"
#include "SVValidation.h"
#include "barcodesExtraction.h"
#include "reverseComplement.h"
#include "help.h"
#include <getopt.h>

// TODO paramètre pour le b de barcodes

int main(int argc, char* argv[]) {
	// Parse arguments
	string bamFile;
	string indexFile;
	string refGenome;
	string outputFile;

	vector<string> s1;
	vector<string> s2;

	unsigned nbBins = 10;
	unsigned maxRegionsLinks = 1000;
	totalCandidates = 0;

	string validCandidatesFile = "candidates.bedpe";
	string cmdLine;


	windowSize = 1000;
	unsigned minVariantSize = 1000;

	mediumVariantsSize = 2000;
	largeVariantsSize = 10000;

	smallVariantsRate = 99;
	mediumVariantsRate = 99;
	largeVariantsRate = 99;

	duplicatesDistance = 10;
	poolSize = 100000;
	minBarcodes = 1;
	unsigned nbThreads = 8;

	const struct option longopts[] = {
		{"bam",					required_argument,	0, 'b'},
		{"index",				required_argument,	0, 'i'},
		{"reference",			required_argument,	0, 'g'},
		{"output",				required_argument,	0, 'o'},

		{"regionSize",			required_argument,	0, 'r'},
		{"minVariantSize",		required_argument,	0, 'v'},

		{"maxLinks",			required_argument,	0, 'n'},
		{"mediumSize",			required_argument,	0, 'M'},
		{"largeSize",			required_argument,	0, 'L'},
		
		{"smallRate",			required_argument,	0, 's'},
		{"mediumRate",			required_argument,	0, 'm'},
		{"largeRate",			required_argument,	0, 'l'},
		
		{"duplicates",			required_argument,	0, 'd'},

		{"threads",				required_argument,	0, 't'},
		{"poolSize",			required_argument,	0, 'p'},
		{"nbBins",				required_argument,	0, 'B'},
		{"minBarcodes",			required_argument,	0, 'c'},
		{"validCandidates",		required_argument,  0,  'C'},
		{0, 0, 0, 0},
	};
	int index;
	int iarg = 0;

	iarg = getopt_long(argc, argv, "b:i:g:r:v:n:M:L:s:m:l:o:d:t:p:B:c:C:", longopts, &index);
	if (iarg == -1) {
		printHelp();
	}
	while (iarg != -1) {
		switch (iarg) {
			case 'b':
				bamFile = optarg;
				break;
			case 'i':
				indexFile = optarg;
				break;
			case 'g':
				refGenome = optarg;
				break;
			case 'r':
				windowSize = stoul(optarg);
				break;
			case 'v':
				minVariantSize = stoul(optarg);
				break;
			case 'n':
				maxRegionsLinks = stoul(optarg);
				break;
			case 'M':
				mediumVariantsSize = stoul(optarg);
				break;
			case 'L':
				largeVariantsSize = stoul(optarg);
				break;
			case 's':
				smallVariantsRate = stod(optarg);
				break;
			case 'm':
				mediumVariantsRate = stod(optarg);
				break;
			case 'l':
				largeVariantsRate = stod(optarg);
				break;
			case 'o':
				outputFile = optarg;
				break;
			case 'd':
				duplicatesDistance = stoul(optarg);
				break;
			case 't':
				nbThreads = stoul(optarg);
				break;
			case 'p':
				poolSize = stoul(optarg);
				break;
			case 'B':
				nbBins = stoul(optarg);
				break;
			case 'c':
				minBarcodes = stoul(optarg);
				break;
			case 'C':
				validCandidatesFile = optarg;
				break;
			default:
				printHelp();
				break;
		}
		iarg = getopt_long(argc, argv, "b:i:g:r:v:n:M:L:s:m:l:o:d:t:p:B:c:C:", longopts, &index);
	}

	if (bamFile.empty() or indexFile.empty() or refGenome.empty() or outputFile.empty()) {
		fprintf(stderr, "Please provide valid input files with options -b, -i, -g and -o.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < argc; i++) {
		cmdLine = cmdLine + argv[i] + " ";
	}

	cerr << "Preparing auxiliary data" << endl;
	prepareAuxiliaryData(bamFile);
	genomeIndex = indexGenome(refGenome);
	
	robin_hood::unordered_map<string*, robin_hood::unordered_map<string*, unsigned>> candidates;
	if (isEmptyFile(validCandidatesFile)) {
		cerr << "Loading the barcodes positions index" << endl;
		BarcodesPositionsIndex barcodesPositionsIndex = loadBarcodesPositionsIndex(indexFile);
		// Remove adjacent positions that correspond to the same window
		removeIndexRedundancy(barcodesPositionsIndex);

		// Process every barcode to look for SV evidence
		cerr << "Computing the number of common barcodes between all the pairs of regions of the genome" << endl;
		totalBarcodes = barcodesPositionsIndex.size();
		candidates = processBarcodes(nbThreads, refIDs, regionsList, nbBins, bamFile, barcodesPositionsIndex, barcodesPositionsIndex.size(), minVariantSize);

		cerr << "Computing and analyzing the distribution of shared barcodes between candidates" << endl;
		Thresholds th = analyzeDistribution(candidates);	
		
		// TODO refaire le message
		cerr << "Removing invalid candidates (regions pairs that share do not share a sufficient number of barcodes and regions that are paired with more than " << maxRegionsLinks << " other regions)" << endl;
		removeInvalidCandidates(candidates, th, maxRegionsLinks);
		// removeCandidatesRegionsLinks(candidates, maxRegionsLinks);
		
		cerr << "Saving all SV candidates to file \"" << validCandidatesFile << "\"" << endl;
		saveSVCandidates(validCandidatesFile, candidates);

		barcodesPositionsIndex.clear();
	} else {
		cerr << "Candidates files \"" << validCandidatesFile << "\" already exists. Loading candidates from it." << endl;
		candidates = loadSVCandidates(validCandidatesFile);
	}


	// Analyze alignments of valid candidates
	cerr << "Number of valid SV candidates to consider : " << totalCandidates << endl;
	robin_hood::unordered_map<pair<string*, string*>, SVSupports, hashPointersPairs> calledSVs = processCandidates(nbThreads, candidates, bamFile);
	candidates.clear();

	robin_hood::unordered_set<StructuralVariant> finalSVs = validatesSVs(calledSVs);

	cerr << "Output " << finalSVs.size() << " SVs" << endl;
	outputSVsAsVCF(finalSVs, cmdLine, refGenome, outputFile);

	// if (remove(validCandidatesFile.c_str()) != 0 ) {
 //    	fprintf(stderr, "Error when deleting file %s.\n", validCandidatesFile.c_str());
 //  	}

	return EXIT_SUCCESS;
}

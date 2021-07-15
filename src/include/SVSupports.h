#ifndef __MYVC_SV_SUPPORTS__
#define __MYVC_SV_SUPPORTS__

/**
	Structure defining the support of a SV
*/
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

// TODO add a function for comparing the support in there

#endif
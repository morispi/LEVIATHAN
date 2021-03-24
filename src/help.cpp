#include "help.h"

void printHelp() {
	printf("%s\n", VERSION);
	printf("Pierre Morisse <pierre.morisse@inria.fr>\n");
	printf("LEVIATHAN: Linked-reads based structural variant caller with barcode indexing\n");
	printf("\n");

	printf("USAGE:\n");
	printf("\tLEVIATHAN -b bamFile.bam -i barcodeIndex.bci -g genome.fasta -o output.vcf [OPTIONS]\n");
	printf("\n");

	printf("INPUT:\n");
	printf("\tbamFile.bam:              BAM file to analyze. Warning: the associated .bai file must exist\n");
	printf("\tbarcodeIndex.bci:         LRez barcode occurrence positions index of the BAM file\n");
	printf("\tgenome.fasta:             The reference genome in FASTA format\n");
	printf("\toutput.vcf:               VCF file where to ouput the SVs\n");
	printf("\n");

	printf("OPTIONS:\n");
	printf("\t-r --regionSize:          Size of the regions on the reference genome to consider (default: 1000)\n");
	printf("\t-v, --minVariantSize:     Minimum size of the SVs to detect (default: same as regionSize)\n");
	printf("\t-n, --maxLinks:           Remove from candidates list all candidates which have a region involved in that much candidates (default: 1000) \n");
	printf("\t-M, --mediumSize:         Minimum size of medium variants (default: 2000)\n");
	printf("\t-L, --largeSize:          Minimum size of large variants (default: 10000)\n");
	printf("\t-s, --smallRate:          Percentile to chose as a threshold in the distribution of the number of shared barcodes for small variants (default: 99)\n");
	printf("\t-m, --mediumRate:         Percentile to chose as a threshold in the distribution of the number of shared barcodes for medium variants (default: 99)\n");
	printf("\t-l, --largeRate:          Percentile to chose as a threshold in the distribution of the number of shared barcodes for large variants (default: 99)\n");
	printf("\t-d, --duplicates:         Consider SV as duplicates if they have the same type and if their breakpoints are within this distance (default: 10)\n");
	printf("\t-t, --threads:            Number of threads (default: 8)\n");
	printf("\t-p, --poolSize:           Size of the thread pool (default: 100000)\n");
	printf("\t-B, --nbBins:             Number of iterations to perform through the barcode index (default: 10)\n");
	printf("\t-c, --minBarcodes:        Always remove candidates that share less than this number of barcodes (default: 1)\n");

	exit(EXIT_FAILURE);
}
# LEVIATHAN

LEVIATHAN (Linked-reads based structural variant caller with barcode indexing) is a structural variant calling method for Linked-Reads data. LEVIATHAN takes as input a BAM file, which can either be generated by a Linked-Reads dedicated mapper such as Long Ranger, or by any other aligner, given the reads are pre-processed in order to extract the barcodes from their sequences and append them to the headers (e.g. using Long Ranger basic). LEVIATHAN output SVs in VCF format.

LEVIATHAN works by, first, indexing the occurrence positions of the barcodes throughout the BAM file using LRez. It then relies on two distinct steps. The first one computes the amount of common barcodes between region pairs of the reference genome, in order to highlight SV candidates, which are pairs of distant regions that share more barcodes than expected. The second step then acts a refining step, and relies on classical short-reads methodologies to further analyze these SV candidates, and determine the types and breakpoints of actual SVs, and filter out false-positives.

Requirements
--------------

  - A Linux based operating system.
  - g++, minimum version 5.5.0.
  - CMake, minimum version 2.8.2.
  - zlib, minimum version 1.2.11.
  
Installation from source
--------------

Clone the LEVIATHAN repository, along with its submodules with:

  ```bash
  git clone --recursive https://github.com/morispi/LEVIATHAN
  ```

Then run the install.sh script:

  ```bash
  ./install.sh
  ```

The installation script will build all dependencies, and build and store the binary in the `bin` folder.

Installation from source
--------------

Alternatively, LEVIATHAN is also distributed as a bioconda package, which can be installed with:

```bash
conda install -c bioconda leviathan
```

Getting started
--------------

An example dataset (BAM and BAM.bai files, as well as reference genome) is provided in the `example` folder.

Please run the following commands to try out LEVIATHAN on this example.

### 1) Build LRez barcode index.

To build the barcode index on the example dataset, run the following command:

`./LRez/bin/LRez index bam -p -b example/example.bam -o example/barcodeIndex.bci`

This should take about 2 seconds, and use 20 KB of RAM.

### 2) Run LEVIATHAN.

To run LEVIATHAN, once the index is buiilt, run the following command:

`./bin/LEVIATHAN -b example/example.bam -i example/barcodeIndex.bci -g example/genome.fasta -o example/SV.vcf`

This should take about 10 seconds, and use at most 100 KB of RAM, using 8 threads.
  
Running LEVIATHAN
--------------

### 1) Build LRez barcode index.

Prior to running LEVIATHAN, the LRez barcode index of the BAM file has to be built. This can be done with the following command:

`./LRez/bin/LRez index bam -p -b bamFile.bam -o barcodeIndex.bci`

  - bamFile.bam:        BAM file to index. Warning: the associated .bai file must exist
  - barcodeIndex.bci:   File where to store the LRez barcode index

### 2) Run LEVIATHAN.

Once the index is built, LEVIATHAN can then be run with the following command:

`./bin/LEVIATHAN -b bamFile.bam -i barcodeIndex.bci -g genome.fasta -o output.vcf`

  - bamFile.bam:        BAM file to analyze. Warning: the associated .bai file must exist
  - barcodeIndex.bci:   LRez barcode index of the BAM file
  - genome.fasta:       Reference genome in FASTA format
  - output.vcf:         VCF file where to ouput the SVs

### Options

      -r --regionSize:          Size of the regions on the reference genome to consider (default: 1000)
      -v, --minVariantSize:     Minimum size of the SVs to detect (default: same as regionSize)
      -n, --maxLinks:           Remove from candidates list all candidates which have a region involved in that much candidates (default: 1000)
      -M, --mediumSize:         Minimum size of medium variants (default: 2000)
      -L, --largeSize:          Minimum size of large variants (default: 10000)
      -s, --smallRate:          Percentile to chose as a threshold in the distribution of the number of shared barcodes for small variants (default: 99)
      -m, --mediumRate:         Percentile to chose as a threshold in the distribution of the number of shared barcodes for medium variants (default: 99)
      -l, --largeRate:          Percentile to chose as a threshold in the distribution of the number of shared barcodes for large variants (default: 99)
      -d, --duplicates:         Consider SV as duplicates if they have the same type and if their breakpoints are within this distance (default: 10)
      -t, --threads:            Number of threads (default: 8)
      -p, --poolSize:           Size of the thread pool (default: 100000)
      -B, --nbBins:             Number of iterations to perform through the barcode index (default: 10)
      -c, --minBarcodes:        Always remove candidates that share less than this number of barcodes (default: 1)
      -C, --candidates:         File where to store valid SV candidates (default: "candidates.bedpe") 

Notes
--------------

LEVIATHAN has been developed and tested on x86-64 GNU/Linux.          
Support for any other platform has not been tested.

Authors
--------------

Pierre Morisse, Fabrice Legeai and Claire Lemaitre.

Reference
--------------

Pierre Morisse, Fabrice Legeai, Claire Lemaitre. LEVIATHAN: efficient discovery of large structural variants by leveraging long-range information from Linked-Reads data. https://doi.org/10.1101/2021.03.25.437002

Contact
--------------

You can report problems and bugs to pierre[dot]morisse[at]inria[dot]fr

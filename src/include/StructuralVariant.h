#ifndef __MYVC_STRUCTURAL_VARIANT__
#define __MYVC_STRUCTURAL_VARIANT__

#include <string>

using namespace std;

struct StructuralVariant {
	string chr1;
	unsigned breakpoint1;
	string chr2;
	unsigned breakpoint2;
	string type;
	unsigned barcodes;
	unsigned pairSupport;

	StructuralVariant() {
		
	}

	StructuralVariant(string c1, unsigned bp1, string c2, unsigned bp2, string t, unsigned b, unsigned p) {
		chr1 = c1;
		breakpoint1 = bp1;
		chr2 = c2;
		breakpoint2 = bp2;
		type = t;
		barcodes = b;
		pairSupport = p;
	}

	bool operator==(const StructuralVariant& other) const{
		return chr1 == other.chr1 and breakpoint1 == other.breakpoint1 and chr2 == other.chr2 and breakpoint2 == other.breakpoint2 and type == other.type and barcodes == other.barcodes and pairSupport == other.pairSupport;
	}
};

namespace std {
	template <> struct hash<StructuralVariant> {
		std::size_t operator()(const StructuralVariant& s) const {
  			using std::size_t;
  			using std::hash;
  			using std::string;
			

  			return hash<string>()(s.chr1) ^ hash<unsigned>()(s.breakpoint1) ^ hash<string>()(s.chr2) ^ hash<unsigned>()(s.breakpoint2) ^ hash<string>()(s.type) ^ hash<unsigned>()(s.barcodes) ^ hash<unsigned>()(s.pairSupport);
		}
	};
}

#endif
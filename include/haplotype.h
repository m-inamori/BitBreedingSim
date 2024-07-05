#ifndef __HAPLOTYPE
#define __HAPLOTYPE

#include <memory>
#include <vector>

class ChromMap;


//////////////////// Haplotype ////////////////////

class Haplotype {
	const std::vector<int>	genos;
	const ChromMap&	chrmap;
	
public:
	Haplotype(const std::vector<int>& gs, const ChromMap& cmap) :
											genos(gs), chrmap(cmap) { }
	
	const std::vector<int>& get_genos() const { return genos; }
	std::size_t size() const { return genos.size(); }
	
public:
	static const Haplotype *reduce(const Haplotype *hap1,
								   const Haplotype *hap2,
								   std::default_random_engine &engine);
	static std::vector<const Haplotype *> create_origin_haplotypes(
											std::size_t num, double len,
											std::size_t num_markers,
											const ChromMap& m);
	static std::vector<int> create_genotypes(std::size_t num_markers);
};

#endif

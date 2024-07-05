#include <random>

#include "../include/haplotype.h"

using namespace std;

const Haplotype *Haplotype::reduce(const Haplotype *hap1,
								   const Haplotype *hap2,
								   std::default_random_engine &engine) {
	// where crossovers happen
	const double	L = 1.0;	// 1 Morgan
	vector<double>	pts;
	double	pt0 = 0.0;
	std::exponential_distribution<>	dist(1.0);
	while(true) {
		const double	pt = pt0 + dist(engine);
		if(pt > L)
			break;
		pts.push_back(pt);
		pt0 = pt;
	}
	
	// 最初はどちらのHaplotypeから取るか
	std::uniform_int_distribution<size_t>	dist_unif(0, 1);
	vector<int>	gts(hap1->size());
	int	ihap = dist_unif(engine);
	size_t	first = 0;
	const Haplotype	*hap = nullptr;
	for(auto p = pts.begin(); p != pts.end(); ++p) {
		hap = ihap == 0 ? hap1 : hap2;
		const size_t	last = static_cast<size_t>(*p * 1000);
		const auto&	geno = hap->get_genos();
		std::copy(geno.begin() + first,
				  geno.begin() + last, gts.begin() + first);
		ihap = ihap == 0 ? 1 : 0;
		first = last;
	}
	hap = ihap == 0 ? hap1 : hap2;
	std::copy(hap->get_genos().begin() + first,
			  hap->get_genos().end(), gts.begin() + first);
	return new Haplotype(gts, hap->chrmap);
}


vector<const Haplotype *> Haplotype::create_origin_haplotypes(
										size_t num, double len,
										size_t num_markers, const ChromMap& m) {
	auto dists = std::make_shared<std::vector<double>>(num_markers);
	for(size_t i = 0; i < num_markers; ++i)
		(*dists)[i] = len * (i + 1) / num_markers;
	
	vector<const Haplotype *>	haplotypes;
	for(size_t i = 0; i < num; ++i) {
		vector<int>	genos = create_genotypes(num_markers);
		haplotypes.push_back(new Haplotype(genos, m));
	}
	return haplotypes;
}

vector<int> Haplotype::create_genotypes(size_t num_markers) {
	std::random_device	seed_gen;
	std::mt19937	engine(seed_gen());
	
	vector<int>	genos(num_markers);
	for(size_t i = 0; i < num_markers; ++i) {
		const std::uint32_t	result = engine();
		genos[i] = (int)(result & 1);
	}
	return genos;
}

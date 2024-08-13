#include <sstream>

#include "../include/BaseInfo.h"
#include "../include/Map.h"
#include "../include/trait.h"
#include "../include/population.h"
#include "../include/common.h"

using namespace std;


//////////////////// BaseInfo ////////////////////

BaseInfo::BaseInfo(const Map *m, const vector<const Trait *>& ts,
											std::uint_fast32_t seed) :
										gmap(m), traits(ts), engine(seed) { }

BaseInfo::~BaseInfo() {
	delete gmap;
	Common::delete_all(traits);
}

const ChromMap& BaseInfo::get_chrom_map(size_t i) const {
	return gmap->get_chr(i);
}

size_t BaseInfo::num_chroms() const {
	return gmap->num_chroms();
}

const string& BaseInfo::get_trait_name(size_t i) const {
	return traits[i]->get_name();
}

BaseInfo *BaseInfo::create_default(int seed) {
	const Map	*gmap = Map::create_default(10, 1000, 1.0, 1000000);
	if(seed == -1) {
		std::random_device	seed_gen;
		return new BaseInfo(gmap, vector<const Trait *>(), seed_gen());
	}
	else {
		return new BaseInfo(gmap, vector<const Trait *>(),
										static_cast<std::uint_fast32_t>(seed));
	}
}

void BaseInfo::set_trait() {
	stringstream	ss;
	ss << "Trait" << num_traits() + 1;
	const Trait	*trait = Trait::create_A_randomly(ss.str(), 1, 0.0,
													1.0, 0.7, gmap, engine);
	traits.push_back(trait);
}

void BaseInfo::set_trait_AD_multi(size_t num_loci, double h2, double H2) {
	stringstream	ss;
	ss << "Trait" << num_traits() + 1;
	const vector<Trait::Locus>	loci = Trait::decide_loci_randomly(
														num_loci, gmap, engine);
	const Trait	*trait = Trait::create_AD_randomly(ss.str(), 0.0, 1.0,
														h2, H2, loci, engine);
	traits.push_back(trait);
}

vector<double> BaseInfo::compute_phenotypes(const Population& pop,
												size_t trait_index) const {
	vector<double>	phenotypes(pop.num_inds());
	const Trait	*trait = traits[trait_index];
	for(size_t i = 0; i < pop.num_inds(); ++i) {
		phenotypes[i] = trait->phenotype(i, pop, engine);
	}
	return phenotypes;
}

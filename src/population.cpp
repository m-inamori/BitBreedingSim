#include <sstream>
#include "../include/population.h"

using namespace std;

//////////////////// ChrPopulation ////////////////////

const ChrPopulation *ChrPopulation::create_origins(size_t num_inds,
													const ChromMap& cmap) {
	std::random_device	seed_gen;
	std::mt19937	engine(seed_gen());
	
	const size_t	num_markers = cmap.num_markers();
	vector<int>	genos(num_markers * num_inds * 2);
	for(size_t i = 0; i < num_markers * num_inds * 2; ++i) {
		const std::uint32_t	result = engine();
		genos[i] = (int)(result & 1);
	}
	return new ChrPopulation(genos, cmap);
}

void ChrPopulation::cross(const vector<Pair>& pairs,
							const ChrPopulation& mathers,
							const ChrPopulation& fathers,
							ChrPopulation& new_population,
							std::mt19937 &engine) {
	for(size_t i = 0; i < new_population.num_inds(); ++i) {
		const size_t	mat_index = pairs[i].first;
		const size_t	pat_index = pairs[i].second;
		mathers.reduce(mat_index, new_population, i, 0, engine);
		fathers.reduce(pat_index, new_population, i, 1, engine);
	}
}

void ChrPopulation::reduce(size_t parent_index,
							ChrPopulation& new_population,
							size_t ind_index, size_t hap_index,
							std::mt19937 &engine) const {
	Iter	new_iter = new_population.get_mut_haplotype(ind_index, hap_index);
	const vector<size_t> pts = chrmap.select_random_crossover_points(engine);
	
	// 最初はどちらのHaplotypeから取るか
	std::uniform_int_distribution<size_t>	dist_unif(0, 1);
	int	ihap = dist_unif(engine);
	size_t	first = 0;
	for(auto p = pts.begin(); p != pts.end(); ++p) {
		ConstIter	iter = get_haplotype(parent_index, ihap);
		const size_t	last = *p;
		std::copy(iter + first, iter + last, new_iter + first);
		ihap = ihap == 0 ? 1 : 0;
		first = last;
	}
	ConstIter	iter = get_haplotype(parent_index, ihap);
	std::copy(iter + first, iter + num_markers(), new_iter + first);
}


//////////////////// Population ////////////////////

Population::~Population() {
	for(auto p = chr_populations.begin(); p != chr_populations.end(); ++p)
		delete *p;
}

const Population *Population::create_origins(size_t num_inds,
									const Map& gmap, const string& name_base) {
	vector<const ChrPopulation *>	chr_pops(gmap.num_chroms());
	for(size_t i = 0; i < gmap.num_chroms(); ++i) {
		chr_pops[i] = ChrPopulation::create_origins(num_inds, *gmap.get_chr(i));
	}
	
	vector<string>	names(num_inds);
	for(size_t j = 0; j < num_inds; ++j) {
		stringstream	ss;
		ss << name_base << j + 1;
		names[j] = ss.str();
	}
	return new Population(chr_pops, gmap, names);
}

Population *Population::cross(size_t num_inds,
						const Population& mothers, const Population& fathers,
						const Map& gmap, const string& name_base) {
	const auto	pairs = make_pairs(num_inds, mothers, fathers);
	
	std::random_device	seed_gen;
	std::mt19937	engine(seed_gen());
	
	// prepare space
	vector<const ChrPopulation *>	chr_pops(gmap.num_chroms());
	for(size_t i = 0; i < gmap.num_chroms(); ++i) {
		const auto&	mat = *mothers.get_chrpops(i);
		const auto&	pat = *fathers.get_chrpops(i);
		auto	*chr_pop = new ChrPopulation(num_inds, *gmap.get_chr(i));
		for(size_t j = 0; j < num_inds; ++j) {
			const size_t	mat_index = pairs[j].first;
			const size_t	pat_index = pairs[j].second;
			mat.reduce(mat_index, *chr_pop, j, 0, engine);
			pat.reduce(pat_index, *chr_pop, j, 1, engine);
		}
		chr_pops[i] = chr_pop;
	}
	
	vector<string>	names(num_inds);
	for(size_t j = 0; j < num_inds; ++j) {
		stringstream	ss;
		ss << name_base << j + 1;
		names[j] = ss.str();
	}
	
	return new Population(chr_pops, gmap, names);
}

vector<Population::Pair> Population::make_pairs(size_t num_inds,
												const Population& mothers,
												const Population& fathers) {
	std::random_device	seed_gen;
	std::mt19937	engine(seed_gen());
	
	std::uniform_int_distribution<size_t>	dist1(0, mothers.num_inds()-1);
	std::uniform_int_distribution<size_t>	dist2(0, fathers.num_inds()-1);
	
	vector<pair<size_t, size_t>>	pairs(num_inds);
	for(size_t i = 0; i < num_inds; ++i) {
		const size_t	mother_index = dist1(engine);
		const size_t	father_index = dist2(engine);
		pairs[i] = make_pair(mother_index, father_index);
	}
	return pairs;
}

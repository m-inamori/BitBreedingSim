#include <sstream>
#include <climits>
#include "../include/population.h"
#include "../include/bitoperation.h"
#include "../include/common.h"

using namespace std;


//////////////////// BitChrPopulation ////////////////////

const BitChrPopulation *BitChrPopulation::create_origins(size_t num_inds,
														const ChromMap& cmap) {
	std::random_device	seed_gen;
	std::mt19937_64	engine(seed_gen());
	std::uniform_int_distribution<Int::ull> dist(0, ULLONG_MAX);
	
	const size_t	num_markers = cmap.num_markers();
	const size_t	num_elements = (num_markers + 63) / 64;
	vector<Int::ull>	genos(num_elements * num_inds * 2);
	// 端数は気にしなくてもよい
	for(size_t i = 0; i < num_elements * num_inds * 2; ++i) {
		genos[i] = dist(engine);
	}
	return new BitChrPopulation(genos, num_inds, cmap);
}

void BitChrPopulation::cross(const vector<Pair>& pairs,
							const BitChrPopulation& mothers,
							const BitChrPopulation& fathers, int T) {
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(i, T, mothers, fathers, pairs, *this);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
			(void *(*)(void *))&cross_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		cross_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
}

void BitChrPopulation::cross_in_thread(void *config) {
	auto	*c = (ConfigThread *)config;
	
	std::random_device	seed_gen;
	std::mt19937	engine(seed_gen());
	
	BitChrPopulation&	pops = c->new_population;
	for(size_t i = c->first; i < pops.get_num_inds(); i += c->num_threads) {
		const size_t	mat_index = c->pairs[i].first;
		const size_t	pat_index = c->pairs[i].second;
		pops.cross_each(c->mothers, c->fathers,
									mat_index, pat_index, i, engine);
	}
}

void BitChrPopulation::cross_each(const BitChrPopulation& mathers,
								  const BitChrPopulation& fathers,
								  size_t mat_index, size_t pat_index,
								  size_t ind_index, std::mt19937 &engine) {
	mathers.reduce(mat_index, *this, ind_index, 0, engine);
	fathers.reduce(pat_index, *this, ind_index, 1, engine);
}

void BitChrPopulation::reduce(size_t parent_index,
							  BitChrPopulation& new_population,
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
		if(last == 0) {
			ihap = ihap == 0 ? 1 : 0;
			continue;
		}
		const size_t	first_q = first / 64;
		const size_t	first_r = first % 64;
		const size_t	last_q = (last - 1) / 64;
		
		const Int::ull	mask = BitOperation::upper_mask(first_r);
		*(new_iter + first_q) &= (~0) ^ mask;
		*(new_iter + first_q) |= (*(iter + first_q)) & mask;
		if(first_q < last_q) {
			std::copy(iter + first_q + 1, iter + last_q + 1,
											new_iter + first_q + 1);
		}
		ihap = ihap == 0 ? 1 : 0;
		first = last;
	}
	ConstIter	iter = get_haplotype(parent_index, ihap);
	const size_t	first_q = first / 64;
	const size_t	first_r = first % 64;
	
	const Int::ull	mask = BitOperation::upper_mask(first_r);
	*(new_iter + first_q) &= (~0) ^ mask;
	*(new_iter + first_q) |= (*(iter + first_q)) & mask;
	std::copy(iter + first_q + 1, iter + num_elements(),
										new_iter + first_q + 1);
}

void BitChrPopulation::write(ostream& os) const {
	// あとでChrの名前もちゃんと入れる
	for(size_t i = 0; i < num_markers(); ++i) {
		os << "Chr1" << '\t' << (i+1)*1000 << '\t' << '.' << '\t'
				<< 'A' << '\t' << 'C' << '\t' << '.' << '\t'
				<< "PASS" << '\t' << '.' << '\t' << "GT";
		for(size_t k = 0; k < num_inds; ++k)
			os << '\t' << get_genotype(k, i);
		os << "\n";
	}
}

string BitChrPopulation::get_genotype(size_t id_ind, size_t id_marker) const {
	size_t	q = id_marker / 64;
	size_t	r = id_marker % 64;
	const size_t	offset1 = id_ind * 2 * num_elements();
	const size_t	offset2 = offset1 + num_elements();
	const Int::ull	gt1 = (genos[offset1+q] >> r) & 1;
	const Int::ull	gt2 = (genos[offset2+q] >> r) & 1;
	stringstream	ss;
	ss << gt1 << '|' << gt2;
	return ss.str();
}


//////////////////// Population ////////////////////

Population::~Population() {
	for(auto p = chr_populations.begin(); p != chr_populations.end(); ++p)
		delete *p;
}

const Population *Population::create_origins(size_t num_inds,
									const Map& gmap, const string& name_base) {
	vector<const BitChrPopulation *>	chr_pops(gmap.num_chroms());
	for(size_t i = 0; i < gmap.num_chroms(); ++i) {
		const auto&	cmap = *gmap.get_chr(i);
		chr_pops[i] = BitChrPopulation::create_origins(num_inds, cmap);
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
						const Map& gmap, const string& name_base, int T) {
	const auto	pairs = make_pairs(num_inds, mothers, fathers);
	
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(i, T, mothers, fathers, pairs);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
			(void *(*)(void *))&cross_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		cross_in_thread(configs[i]);
#endif
	
	// make them to be const
	const auto&	cpops = configs[0]->chr_pops;
	const vector<const BitChrPopulation *>	chr_pops(cpops.begin(),
														cpops.end());
	
	Common::delete_all(configs);
	
	vector<string>	names(num_inds);
	for(size_t j = 0; j < num_inds; ++j) {
		stringstream	ss;
		ss << name_base << j + 1;
		names[j] = ss.str();
	}
	
	return new Population(chr_pops, gmap, names);
}

#include <iostream>

void Population::cross_in_thread(void *config) {
	auto	*c = (ConfigThread *)config;
	
	for(size_t i = c->first; i < c->num_chroms(); i += c->num_threads) {
		const auto&	mat = *c->mothers.get_chrpops(i);
		const auto&	pat = *c->fathers.get_chrpops(i);
		const ChromMap&	cmap = c->mothers.get_chrmap(i);
		auto	*chr_pop = new BitChrPopulation(c->num_inds(), cmap);
		chr_pop->cross(c->pairs, mat, pat, 1);
		c->chr_pops[i] = chr_pop;
	}
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

void Population::write(ostream& os) const {
	write_header(os);
	for(auto p = chr_populations.begin(); p != chr_populations.end(); ++p)
		(*p)->write(os);
}

void Population::write_header(ostream& os) const {
	os << "##fileformat=VCFv4.2\n";
	os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	os << "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT";
	for(auto p = names.begin(); p != names.end(); ++p) {
		os << '\t' << *p;
	}
	os << "\n";
}

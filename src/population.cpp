#include <sstream>
#include <climits>
#include "../include/BaseInfo.h"
#include "../include/trait.h"
#include "../include/population.h"
#include "../include/VCF.h"
#include "../include/bitoperation.h"
#include "../include/common.h"

using namespace std;


//////////////////// BitChrPopulation ////////////////////

const BitChrPopulation *BitChrPopulation::create_origins(size_t num_inds,
													const ChromMap& cmap,
													std::mt19937_64& engine) {
	std::uniform_int_distribution<Int::ull> dist(0, ULLONG_MAX);
	
	const size_t	num_markers = cmap.get_num_markers();
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
								const BitChrPopulation& fathers,
								std::uint_fast32_t seed0) {
	std::mt19937	engine(seed0);
	for(size_t ind_index = 0; ind_index < num_inds; ++ind_index) {
		const size_t	mat_index = pairs[ind_index].first;
		const size_t	pat_index = pairs[ind_index].second;
		mothers.reduce(mat_index, *this, ind_index, 0, engine);
		fathers.reduce(pat_index, *this, ind_index, 1, engine);
	}
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
	for(size_t i = 0; i < num_markers(); ++i) {
		vector<string>	gts(num_inds);
		for(size_t k = 0; k < num_inds; ++k)
			gts[k] = get_genotype(k, i);
		VCF::write_data_line(os, chrmap.get_name(), 
								chrmap.get_position(i), gts);
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

int BitChrPopulation::get_int_genotype(size_t id_ind, size_t id_marker) const {
	size_t	q = id_marker / 64;
	size_t	r = id_marker % 64;
	const size_t	offset1 = id_ind * 2 * num_elements();
	const size_t	offset2 = offset1 + num_elements();
	const Int::ull	gt1 = (genos[offset1+q] >> r) & 1;
	const Int::ull	gt2 = (genos[offset2+q] >> r) & 1;
	return static_cast<int>(gt1 + gt2) - 1;
}

BitChrPopulation *BitChrPopulation::select(
							const vector<size_t>& indices) const {
	const size_t	N = num_elements() * 2;		// data per individual
	vector<Int::ull>	selected_genos(N * indices.size());
	for(size_t i = 0; i < indices.size(); ++i) {
		const size_t	ind_id = indices[i];
		std::copy(genos.begin() + N*ind_id, genos.begin() + N*(ind_id+1),
											selected_genos.begin() + N*i);
	}
	return new BitChrPopulation(selected_genos, indices.size(), chrmap);
}


//////////////////// Population ////////////////////

Population::~Population() {
	for(auto p = chr_populations.begin(); p != chr_populations.end(); ++p)
		delete *p;
}

void Population::set_phenotypes(const BaseInfo *info) {
	for(size_t i = 0; i < info->num_traits(); ++i) {
		const auto	pheno = info->compute_phenotypes(*this, i);
		phenotypes.push_back(pheno);
		traits.push_back(info->get_trait(i));
	}
}

vector<vector<double>> Population::compute_phenotypes(const BaseInfo *info,
														std::mt19937 &engine) {
	vector<vector<double>>	phenotypes(info->num_traits());
	for(size_t i = 0; i < info->num_traits(); ++i) {
		const Trait	*trait = info->get_trait(i);
		for(size_t k = 0; k < num_inds(); ++k) {
			phenotypes[i][k] = trait->phenotype(k, *this, engine);
		}
	}
	return phenotypes;
}

const Population *Population::create_origins(size_t num_inds,
							const BaseInfo *info, const string& name_base) {
	const Map&	gmap = info->get_map();
	std::mt19937&	engine = info->get_random_engine();
	std::mt19937_64	engine64(engine());
	vector<const BitChrPopulation *>	chr_pops(gmap.num_chroms());
	for(size_t i = 0; i < gmap.num_chroms(); ++i) {
		const auto&	cmap = gmap.get_chr(i);
		chr_pops[i] = BitChrPopulation::create_origins(num_inds,
														cmap, engine64);
	}
	
	vector<string>	names(num_inds);
	for(size_t j = 0; j < num_inds; ++j) {
		stringstream	ss;
		ss << name_base << j + 1;
		names[j] = ss.str();
	}
	
	Population	*pop = new Population(chr_pops, gmap, names);
	pop->set_phenotypes(info);
	return pop;
}

Population *Population::cross(size_t num_inds,
						const Population& mothers, const Population& fathers,
						const BaseInfo *info, const string& name_base, int T) {
	std::mt19937&	engine = info->get_random_engine();
	const auto	pairs = make_pairs(num_inds, mothers, fathers, engine);
for(const auto& p : pairs) cout << p.first << ',' << p.second << endl;
	const std::uint_fast32_t	seed = engine();
	vector<BitChrPopulation *>	chr_pops_(info->num_chroms());
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(i, T, mothers, fathers,
												pairs, chr_pops_, seed);
	
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
	const vector<BitChrPopulation *>	chr_ps = configs[0]->chr_pops;
	const vector<const BitChrPopulation *>	chr_pops(chr_ps.begin(),
															chr_ps.end());
	
	Common::delete_all(configs);
	
	vector<string>	names(num_inds);
	for(size_t j = 0; j < num_inds; ++j) {
		stringstream	ss;
		ss << name_base << j + 1;
		names[j] = ss.str();
	}
	
	const Map&	gmap = info->get_map();
	Population	*pop = new Population(chr_pops, gmap, names);
	pop->set_phenotypes(info);
	return pop;
}

void Population::cross_in_thread(void *config) {
	auto	*c = (ConfigThread *)config;
	
	for(size_t i = c->first; i < c->num_chroms(); i += c->num_threads) {
		const auto&	mat = *c->mothers.get_chrpops(i);
		const auto&	pat = *c->fathers.get_chrpops(i);
		const ChromMap&	cmap = c->mothers.get_chrmap(i);
		auto	*chr_pop = new BitChrPopulation(c->num_inds(), cmap);
		chr_pop->cross(c->pairs, mat, pat, c->seed0 + i);
		c->chr_pops[i] = chr_pop;
	}
}

vector<Population::Pair> Population::make_pairs(size_t num_inds,
												const Population& mothers,
												const Population& fathers,
												std::mt19937& engine) {
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
	VCF::write_header(os, names);
	for(auto p = chr_populations.begin(); p != chr_populations.end(); ++p)
		(*p)->write(os);
}

void Population::write_phenotypes(ostream& os) const {
	os << "name";
	for(const auto& trait : traits) {
		os << ',' << trait->get_name();
	}
	os << "\n";
	for(size_t i = 0; i < num_inds(); ++i) {
		os << names[i];
		for(size_t k = 0; k < traits.size(); ++k) {
			os << ',' << phenotypes[k][i];
		}
		os << "\n";
	}
}

double Population::mean(size_t i) const {
	return std::accumulate(phenotypes[i].begin(),
							phenotypes[i].end(), 0.0) / num_inds();
}

double Population::stddev(std::size_t i) const {
	const double	m = mean(i);
	double	s2 = 0;
	for(const auto& p : phenotypes[i]) {
		s2 += (p - m) * (p - m);
	}
	return sqrt(s2 / num_inds());
}

void Population::dispay_QTLs(size_t i) const {
	const Trait	*trait = traits[i];
	const vector<Trait::Locus>	loci = trait->get_loci();
	const vector<double>	additives = trait->get_addivtives();
	const vector<double>	dominants = trait->get_dominants();
	
	cout << "Chr";
	for(size_t i = 0; i < trait->num_QTLs(); ++i)
		cout << ',' << loci[i].first + 1;
	cout << endl;
	cout << "pos";
	for(size_t i = 0; i < trait->num_QTLs(); ++i)
		cout << ',' << loci[i].second + 1;
	cout << endl;
	cout << "additive";
	for(size_t i = 0; i < trait->num_QTLs(); ++i)
		cout << ',' << additives[i];
	cout << endl;
	cout << "dominant";
	for(size_t i = 0; i < trait->num_QTLs(); ++i)
		cout << ',' << dominants[i];
	cout << endl;
	
	for(size_t k = 0; k < num_inds(); ++k) {
		cout << names[k];
		for(size_t i = 0; i < trait->num_QTLs(); ++i)
			cout << ',' << get_genotype(k, loci[i].first, loci[i].second);
		cout << endl;
	}
}

Population *Population::select(const vector<size_t>& indices) const {
	vector<const BitChrPopulation *>	selected_chr_pops(num_chroms());
	for(size_t i = 0; i < num_chroms(); ++i) {
		selected_chr_pops[i] = chr_populations[i]->select(indices);
	}
	
	vector<string>	selected_names(indices.size());
	for(size_t i = 0; i < indices.size(); ++i) {
		selected_names[i] = names[indices[i]];
	}
	Population	*selected_pop = new Population(selected_chr_pops,
														gmap, selected_names);
	
	selected_pop->phenotypes.resize(phenotypes.size());
	for(size_t i = 0; i < phenotypes.size(); ++i) {
		selected_pop->phenotypes[i] = select_phenotypes(indices, i);
	}
	selected_pop->traits = traits;
	
	return selected_pop;
}

vector<double> Population::select_phenotypes(const vector<size_t>& indices,
															size_t i) const {
	vector<double>	selected_phenotypes(indices.size());
	for(size_t k = 0; k < indices.size(); ++k) {
		selected_phenotypes[k] = phenotypes[i][indices[k]];
	}
	return selected_phenotypes;
}

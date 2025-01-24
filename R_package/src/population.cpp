#include <fstream>
#include <sstream>
#include <climits>
#include <Rcpp.h>
#include "BaseInfo.h"
#include "trait.h"
#include "population.h"
#include "VCF.h"
#include "bitarray.h"
#include "bitoperation.h"
#include "common.h"

using namespace std;
using namespace Rcpp;


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
		BitArray::copy(iter, first, last, new_iter);
		ihap = ihap == 0 ? 1 : 0;
		first = last;
	}
	ConstIter	iter = get_haplotype(parent_index, ihap);
	BitArray::copy(iter, first, num_elements() * 64, new_iter);
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
	ConstIter	iter1 = get_haplotype(id_ind, 0);
	ConstIter	iter2 = get_haplotype(id_ind, 1);
	const int	gt1 = BitArray::get(iter1, id_marker);
	const int	gt2 = BitArray::get(iter2, id_marker);
	stringstream	ss;
	ss << gt1 << '|' << gt2;
	return ss.str();
}

int BitChrPopulation::get_int_genotype(size_t id_ind, size_t id_marker) const {
	ConstIter	iter1 = get_haplotype(id_ind, 0);
	ConstIter	iter2 = get_haplotype(id_ind, 1);
	const int	gt1 = BitArray::get(iter1, id_marker);
	const int	gt2 = BitArray::get(iter2, id_marker);
	return gt1 + gt2 - 1;
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

const BitChrPopulation *BitChrPopulation::join(const BitChrPopulation *pop1,
											   const BitChrPopulation *pop2) {
	vector<Int::ull>	genos;
	Common::connect_vector(pop1->genos, pop2->genos, genos);
	const size_t	num_inds = pop1->num_inds + pop2->num_inds;
	return new BitChrPopulation(genos, num_inds, pop1->chrmap);
}


//////////////////// Population ////////////////////

Population::~Population() {
	for(auto p = chr_populations.begin(); p != chr_populations.end(); ++p)
		delete *p;
}

size_t Population::num_markers() const {
	size_t	num = 0;
	for(const auto& chr_pop : chr_populations) {
		num += chr_pop->num_markers();
	}
	return num;
}

size_t Population::num_traits() const {
	return info->num_traits();
}

const ChromMap&	Population::get_chrmap(std::size_t i) const {
	return info->get_chrom_map(i);
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
	
	vector<string>	mats(num_inds, "0");
	vector<string>	pats(num_inds, "0");
	
	Population	*pop = new Population(chr_pops, info, names, mats, pats);
	pop->set_phenotypes(info);
	return pop;
}

vector<Population::Pair> Population::make_pairs_randomly(size_t num_inds,
												const Population& mothers,
												const Population& fathers,
												std::mt19937& engine) {
	std::uniform_int_distribution<size_t>	dist1(0, mothers.num_inds()-1);
	std::uniform_int_distribution<size_t>	dist2(0, fathers.num_inds()-1);
	
	vector<Pair>	pairs(num_inds);
	for(size_t i = 0; i < num_inds; ++i) {
		const size_t	mother_index = dist1(engine);
		const size_t	father_index = dist2(engine);
		pairs[i] = make_pair(mother_index, father_index);
	}
	return pairs;
}

vector<Population::Pair> Population::make_pairs_by_table(
												const vector<Triplet>& table,
												const Population& mothers,
												const Population& fathers) {
	map<string, size_t>	dic_mat;
	for(size_t i = 0; i < mothers.names.size(); ++i) {
		dic_mat.insert(make_pair(mothers.names[i], i));
	}
	map<string, size_t>	dic_pat;
	for(size_t i = 0; i < fathers.names.size(); ++i) {
		dic_pat.insert(make_pair(fathers.names[i], i));
	}
	
	size_t	num_inds = 0;
	for(const auto& t : table) {
		num_inds += get<2>(t);
	}
	vector<Pair>	pairs(num_inds);
	size_t	i = 0;
	for(const auto& t: table) {
		const string&	mat = get<0>(t);
		const string&	pat = get<1>(t);
		const size_t&	num = get<2>(t);
		auto	p = dic_mat.find(mat);
		auto	q = dic_pat.find(pat);
		if(p == dic_mat.end() || q == dic_pat.end()) {
			// If an error is detected, it is returned to R immediately.
			// Error details are generated and output on the R side.
			stop("parent in table not in parent population.");
		}
		for(size_t k = 0; k < num; ++k) {
			pairs[i] = make_pair(p->second, q->second);
			++i;
		}
	}
	return pairs;
}

Population *Population::cross_randomly(size_t num_inds,
						const Population& mothers, const Population& fathers,
						const BaseInfo *info, const string& name_base, int T) {
	std::mt19937&	engine = info->get_random_engine();
	const auto	pairs = make_pairs_randomly(num_inds, mothers, fathers, engine);
	return cross(pairs, mothers, fathers, info, name_base, engine, T);
}

Population *Population::cross_by_table(const vector<Triplet>& table,
						const Population& mothers, const Population& fathers,
						const BaseInfo *info, const string& name_base, int T) {
	std::mt19937&	engine = info->get_random_engine();
	const auto	pairs = make_pairs_by_table(table, mothers, fathers);
	return cross(pairs, mothers, fathers, info, name_base, engine, T);
}

Population *Population::cross(const vector<Pair>& pairs,
						const Population& mothers, const Population& fathers,
						const BaseInfo *info, const string& name_base,
						std::mt19937& engine, int T) {
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
	
	const size_t	num_inds = pairs.size();
	vector<string>	names(num_inds);
	for(size_t j = 0; j < num_inds; ++j) {
		stringstream	ss;
		ss << name_base << j + 1;
		names[j] = ss.str();
	}
	
	vector<string>	mats(num_inds);
	vector<string>	pats(num_inds);
	std::transform(pairs.begin(), pairs.end(), mats.begin(),
			[&mothers](const Pair& e) { return mothers.get_name(e.first); });
	std::transform(pairs.begin(), pairs.end(), pats.begin(),
			[&fathers](const Pair& e) { return fathers.get_name(e.second); });
	
	Population	*pop = new Population(chr_pops, info, names, mats, pats);
	pop->set_phenotypes(info);
	return pop;
}

void Population::cross_in_thread(void *config) {
	auto	*c = (ConfigThread *)config;
	
	for(size_t i = c->first; i < c->num_chroms(); i += c->num_threads) {
		const auto&	mat = *c->mothers.get_chrpop(i);
		const auto&	pat = *c->fathers.get_chrpop(i);
		const ChromMap&	cmap = c->mothers.get_chrmap(i);
		auto	*chr_pop = new BitChrPopulation(c->num_inds(), cmap);
		chr_pop->cross(c->pairs, mat, pat, c->seed0 + i);
		c->chr_pops[i] = chr_pop;
	}
}

void Population::write_VCF(ostream& os) const {
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
	vector<const BitChrPopulation *>	new_chr_pops(num_chroms());
	for(size_t i = 0; i < num_chroms(); ++i) {
		new_chr_pops[i] = chr_populations[i]->select(indices);
	}
	
	vector<string>	new_names(indices.size());
	for(size_t i = 0; i < indices.size(); ++i) {
		new_names[i] = names[indices[i]];
	}
	vector<string>	new_mats(indices.size());
	vector<string>	new_pats(indices.size());
	std::transform(indices.begin(), indices.end(), new_mats.begin(),
							[this](const size_t& i) { return this->mats[i]; });
	std::transform(indices.begin(), indices.end(), new_mats.begin(),
							[this](const size_t& i) { return this->pats[i]; });
	
	Population	*new_pop = new Population(new_chr_pops, info,
												new_names, new_mats, new_pats);
	
	new_pop->phenotypes.resize(phenotypes.size());
	for(size_t i = 0; i < phenotypes.size(); ++i) {
		new_pop->phenotypes[i] = select_phenotypes(indices, i);
	}
	new_pop->traits = traits;
	
	return new_pop;
}

const Population *Population::join(const Population *pop1,
								   const Population *pop2) {
	vector<const BitChrPopulation *>	chr_pops(pop1->num_chroms());
	for(size_t i = 0; i < pop1->num_chroms(); ++i) {
		chr_pops[i] = BitChrPopulation::join(pop1->chr_populations[i],
											 pop2->chr_populations[i]);
	}
	vector<string>	names;
	Common::connect_vector(pop1->names, pop2->names, names);
	vector<string>	mats;
	Common::connect_vector(pop1->mats, pop2->mats, mats);
	vector<string>	pats;
	Common::connect_vector(pop1->pats, pop2->pats, pats);
	auto	*new_pop = new Population(chr_pops, pop1->info, names, mats, pats);
	
	const size_t	N = pop1->traits.size();
	vector<vector<double>>	phenos(N);
	for(size_t i = 0; i < N; ++i) {
		Common::connect_vector(pop1->phenotypes[i],
							   pop2->phenotypes[i], phenos[i]);
	}
	new_pop->phenotypes = phenos;
	new_pop->traits = pop1->traits;
	return new_pop;
}

vector<double> Population::select_phenotypes(const vector<size_t>& indices,
															size_t i) const {
	vector<double>	selected_phenotypes(indices.size());
	for(size_t k = 0; k < indices.size(); ++k) {
		selected_phenotypes[k] = phenotypes[i][indices[k]];
	}
	return selected_phenotypes;
}


// [[Rcpp::export]]
SEXP createOrigins(SEXP num_inds, SEXP info, SEXP name_base) {
	size_t num_inds_cpp = as<size_t>(num_inds);
	Rcpp::XPtr<BaseInfo> info_cpp(info);
	std::string name_base_cpp = as<std::string>(name_base);
	Rcpp::XPtr<Population> ptr(const_cast<Population*>(
						Population::create_origins(
						num_inds_cpp, info_cpp.get(), name_base_cpp)), true);
	ptr.attr("class") = "Population";
	return ptr;
}

// [[Rcpp::export]]
SEXP crossPopsRandomly(SEXP num_inds, SEXP mothers, SEXP fathers,
						SEXP name_base, Rcpp::CharacterVector names, int T) {
	// あとで使う
	vector<string>	names_cpp(names.size());
	for(int i = 0; i < names.size(); ++i) {
		names_cpp[i] = Rcpp::as<string>(names[i]);
	}
	
	size_t num_inds_cpp = as<size_t>(num_inds);
	Rcpp::XPtr<Population> mothers_cpp(mothers);
	Rcpp::XPtr<Population> fathers_cpp(fathers);
	const BaseInfo	*info = mothers_cpp.get()->get_info();
	const std::string name_base_cpp = as<std::string>(name_base);
	Rcpp::XPtr<Population> ptr(const_cast<Population *>(
						Population::cross_randomly(num_inds_cpp,
										*mothers_cpp.get(),
										*fathers_cpp.get(), info,
										name_base_cpp, T)), true);
	ptr.attr("class") = "Population";
	return ptr;
}

// [[Rcpp::export]]
SEXP crossPopsByTable(DataFrame df, SEXP mothers, SEXP fathers,
											SEXP name_base, int T) {
	CharacterVector mat = df["mats"];
	CharacterVector pat = df["pats"];
	IntegerVector num = df["nums"];
	
	int n = df.nrows();
	
	vector<Population::Triplet>	table;
	for (int i = 0; i < n; ++i) {
		std::string mat_str = as<std::string>(mat[i]);
		std::string pat_str = as<std::string>(pat[i]);
		std::size_t num_val = static_cast<std::size_t>(num[i]);
		
		table.push_back(std::make_tuple(mat_str, pat_str, num_val));
	}
	
	Rcpp::XPtr<Population> mothers_cpp(mothers);
	Rcpp::XPtr<Population> fathers_cpp(fathers);
	const BaseInfo	*info = mothers_cpp.get()->get_info();
	std::string name_base_cpp = as<std::string>(name_base);
	Rcpp::XPtr<Population> ptr(const_cast<Population*>(
						Population::cross_by_table(table,
											*mothers_cpp.get(),
											*fathers_cpp.get(), info,
											name_base_cpp, T)), true);
	ptr.attr("class") = "Population";
	return ptr;
}

// [[Rcpp::export]]
void writeVCF(SEXP pop, const std::string& filename) {
	Rcpp::XPtr<Population> popCpp(pop);
	std::ofstream	ofs(filename);
	if(!ofs.is_open()) {
		Rcpp::stop("Failed to open file: " + filename);
	}
	popCpp.get()->write_VCF(ofs);
	ofs.close();
}

// [[Rcpp::export]]
int getNumInds(SEXP pop) {
	Rcpp::XPtr<Population> pop_cpp(pop);
	return static_cast<int>(pop_cpp.get()->num_inds());
}

// [[Rcpp::export]]
int getNumChromsPop(SEXP pop) {
	Rcpp::XPtr<Population> pop_cpp(pop);
	return static_cast<int>(pop_cpp.get()->num_chroms());
}

// [[Rcpp::export]]
SEXP getPhenotypesCpp(SEXP pop, int i) {
	Rcpp::XPtr<Population> pop_cpp(pop);
	const auto	phenos = pop_cpp.get()->get_phenotypes(i-1);
	NumericVector rphenos(phenos.begin(), phenos.end());
	return rphenos;
}

// [[Rcpp::export]]
SEXP selectPop(SEXP pop, NumericVector indices_R) {
	Rcpp::XPtr<Population> pop_cpp(pop);
	vector<size_t> indices_cpp(indices_R.size());
    for(int i = 0; i < indices_R.size(); ++i) {
        indices_cpp[i] = static_cast<size_t>(indices_R[i]-1);
    }
	Rcpp::XPtr<Population> ptr(const_cast<Population*>(
								pop_cpp.get()->select(indices_cpp)), true);
	return ptr;
}

// [[Rcpp::export]]
NumericMatrix getGenotypes(SEXP pop) {
	Rcpp::XPtr<Population> pop_cpp(pop);
	const Population	*ptr_pop = pop_cpp.get();
	NumericMatrix mat(ptr_pop->num_inds(), ptr_pop->num_markers());
	size_t	mat_index = 0;
	for(size_t i = 0; i < ptr_pop->num_chroms(); ++i) {
		const auto	*chr_pop = ptr_pop->get_chrpop(i);
		const vector<Int::ull>&	genos = chr_pop->get_genos();
		for(size_t ind_id = 0; ind_id < chr_pop->get_num_inds(); ++ind_id) {
			const size_t	num_elems = chr_pop->num_elements();
			for(size_t gid = 0; gid < num_elems; ++gid) {
				const Int::ull&	geno1 = genos[num_elems*2*ind_id+gid];
				const Int::ull&	geno2 = genos[num_elems*(2*ind_id+1)+gid];
				const size_t	num = gid != num_elems - 1 ? 64 :
											chr_pop->num_markers() - 64 * gid;
				for(size_t k = 0; k < num; ++k) {
					const int	gt = static_cast<int>(((geno1 >> k) & 1) +
													  ((geno2 >> k) & 1)) - 1;
					mat(ind_id, mat_index + (gid<<6) + k) = gt;
				}
			}
		}
		mat_index += chr_pop->num_markers();
	}
	
	CharacterVector colnames(ptr_pop->num_markers());
	for(size_t i = 0; i < ptr_pop->num_markers(); ++i) {
		stringstream	ss;
		ss << "marker" << std::setw(8) << std::setfill('0') << i + 1;
		colnames[i] = ss.str();
	}
	const vector<string>&	names = ptr_pop->get_names();
	CharacterVector rownames(names.begin(), names.end());
	mat.attr("dimnames") = List::create(rownames, colnames);
	return mat;
}

// [[Rcpp::export]]
NumericMatrix getGenotypes_naive(SEXP pop) {
	Rcpp::XPtr<Population> pop_cpp(pop);
	const Population	*ptr_pop = pop_cpp.get();
	NumericMatrix mat(ptr_pop->num_inds(), ptr_pop->num_markers());
	size_t	total_marker_index = 0;
	for(size_t i = 0; i < ptr_pop->num_chroms(); ++i) {
		const auto	*chrpop = ptr_pop->get_chrpop(i);
		for(size_t m_id = 0; m_id < chrpop->num_markers(); ++m_id) {
			for(size_t ind_id = 0; ind_id < chrpop->get_num_inds(); ++ind_id) {
				const int	gt = ptr_pop->get_int_genotype(ind_id, i, m_id);
				mat(ind_id, total_marker_index) = gt;
			}
			total_marker_index += 1;
		}
	}
	return mat;
}

// [[Rcpp::export]]
CharacterMatrix getPhasedGenotypes(SEXP pop) {
	Rcpp::XPtr<Population> pop_cpp(pop);
	const Population	*ptr_pop = pop_cpp.get();
	CharacterMatrix mat(ptr_pop->num_inds(), ptr_pop->num_markers());
	size_t	mat_index = 0;
	for(size_t i = 0; i < ptr_pop->num_chroms(); ++i) {
		const auto	*chr_pop = ptr_pop->get_chrpop(i);
		const vector<Int::ull>&	genos = chr_pop->get_genos();
		for(size_t ind_id = 0; ind_id < chr_pop->get_num_inds(); ++ind_id) {
			const size_t	num_elems = chr_pop->num_elements();
			for(size_t gid = 0; gid < num_elems; ++gid) {
				const Int::ull&	geno1 = genos[num_elems*2*ind_id+gid];
				const Int::ull&	geno2 = genos[num_elems*(2*ind_id+1)+gid];
				const size_t	num = gid != num_elems - 1 ? 64 :
											chr_pop->num_markers() - 64 * gid;
				for(size_t k = 0; k < num; ++k) {
					stringstream	ss;
					ss << ((geno1 >> k) & 1) << '|' << ((geno2 >> k) & 1);
					mat(ind_id, mat_index + (gid<<6) + k) = ss.str();
				}
			}
		}
		mat_index += chr_pop->num_markers();
	}
	
	CharacterVector colnames(ptr_pop->num_markers());
	for(size_t i = 0; i < ptr_pop->num_markers(); ++i) {
		stringstream	ss;
		ss << "marker" << std::setw(8) << std::setfill('0') << i + 1;
		colnames[i] = ss.str();
	}
	const vector<string>&	names = ptr_pop->get_names();
	CharacterVector rownames(names.begin(), names.end());
	mat.attr("dimnames") = List::create(rownames, colnames);
	return mat;
}

// [[Rcpp::export]]
NumericMatrix getPhasedIntGenotypes(SEXP pop) {
	Rcpp::XPtr<Population> pop_cpp(pop);
	const Population	*ptr_pop = pop_cpp.get();
	NumericMatrix mat(ptr_pop->num_inds()*2, ptr_pop->num_markers());
	size_t	mat_index = 0;
	for(size_t i = 0; i < ptr_pop->num_chroms(); ++i) {
		const auto	*chr_pop = ptr_pop->get_chrpop(i);
		const vector<Int::ull>&	genos = chr_pop->get_genos();
		for(size_t ind_id = 0; ind_id < chr_pop->get_num_inds(); ++ind_id) {
			const size_t	num_elems = chr_pop->num_elements();
			for(size_t gid = 0; gid < num_elems; ++gid) {
				const Int::ull&	geno1 = genos[num_elems*2*ind_id+gid];
				const Int::ull&	geno2 = genos[num_elems*(2*ind_id+1)+gid];
				const size_t	num = gid != num_elems - 1 ? 64 :
											chr_pop->num_markers() - 64 * gid;
				for(size_t k = 0; k < num; ++k) {
					mat(ind_id*2,   mat_index + (gid<<6) + k) = (geno1>>k)&1;
					mat(ind_id*2+1, mat_index + (gid<<6) + k) = (geno2>>k)&1;
				}
			}
		}
		mat_index += chr_pop->num_markers();
	}
	
	CharacterVector colnames(ptr_pop->num_markers());
	for(size_t i = 0; i < ptr_pop->num_markers(); ++i) {
		stringstream	ss;
		ss << "marker" << std::setw(8) << std::setfill('0') << i + 1;
		colnames[i] = ss.str();
	}
	const vector<string>&	names = ptr_pop->get_names();
	CharacterVector rownames(names.size()*2);
	for(size_t i = 0; i < names.size(); ++i) {
		rownames[i*2]   = names[i] + "_mat";
		rownames[i*2+1] = names[i] + "_pat";
	}
	mat.attr("dimnames") = List::create(rownames, colnames);
	return mat;
}

// [[Rcpp::export]]
SEXP getPopulationInfo(SEXP pop) {
	Rcpp::XPtr<Population>	ptr_pop(pop);
	
	Rcpp::List	pop_list = Rcpp::List::create(
		_["num_inds"] = ptr_pop->num_inds(),
		_["num_chroms"] = ptr_pop->num_chroms(),
		_["num_markers"] = ptr_pop->num_markers(),
		_["num_traits"] = ptr_pop->num_traits()
	);
	return pop_list;
}

// [[Rcpp::export]]
Rcpp::DataFrame createNameDataFromPop(SEXP pop) {
	Rcpp::XPtr<Population>	ptr_pop(pop);
	const auto	*p = ptr_pop.get();
    return Rcpp::DataFrame::create(Rcpp::Named("name") = p->get_names(),
                                   Rcpp::Named("mat") = p->get_mats(),
                                   Rcpp::Named("pat") = p->get_pats());
}

// [[Rcpp::export]]
SEXP joinPop(SEXP pop1, SEXP pop2) {
	Rcpp::XPtr<Population> pop1_cpp(pop1);
	Rcpp::XPtr<Population> pop2_cpp(pop2);
	Rcpp::XPtr<Population> ptr(const_cast<Population *>(
							Population::join(pop1_cpp.get(), pop2_cpp.get())));
	return ptr;
}

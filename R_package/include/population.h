#ifndef __POPULATION
#define __POPULATION

#include <ostream>
#include <vector>
#include <random>
#include <cstdint>
#include <Rcpp.h>
#include "Map.h"
#include "GenomicsCommon.h"

class BaseInfo;
class Trait;
class VCF;

namespace GC = GenomicsCommon;


//////////////////// BitChrPopulation ////////////////////

class BitChrPopulation {
	using ConstIter = std::vector<uint64_t>::const_iterator;
	using Iter = std::vector<uint64_t>::iterator;
	using Pair = std::pair<std::size_t, std::size_t>;
	
private:
	std::vector<uint64_t>		genos;
	const std::vector<GC::Pos>&	positions;
	const std::size_t			num_inds;
	const ChromMap&				chrmap;
	
public:
	BitChrPopulation(const std::vector<uint64_t>& gs,
									const std::vector<GC::Pos>& ps,
									std::size_t n,
									const ChromMap& cmap) :
					genos(gs), positions(ps), num_inds(n), chrmap(cmap) { }
	
	BitChrPopulation(std::size_t num_inds_, const std::vector<GC::Pos>& ps,
														const ChromMap& cmap) :
							genos((ps.size()+63)/64*2*num_inds_),
							positions(ps), num_inds(num_inds_), chrmap(cmap) { }
	
	std::size_t num_markers() const { return positions.size(); }
	std::size_t num_elements() const { return (num_markers()+63)/64; }
	std::size_t get_num_inds() const { return num_inds; }
	std::vector<uint64_t> get_genos() const { return genos; }
	std::string get_genotype(std::size_t id_ind, std::size_t id_marker) const;
	int get_int_genotype(std::size_t id_ind, std::size_t id_marker) const;
	int get_haplotype(std::size_t id_ind,
						std::size_t id_marker, std::size_t parent) const;
	ConstIter get_haplotype(std::size_t ind_index, std::size_t hap_id) const {
		return genos.begin() + num_elements() * (ind_index * 2 + hap_id);
	}
	Iter get_mut_haplotype(std::size_t ind_index, std::size_t hap_id) {
		return genos.begin() + num_elements() * (ind_index * 2 + hap_id);
	}
	const std::vector<GC::Pos>& get_positions() const { return positions; }
	
	void cross(const std::vector<Pair>& pairs,
				const BitChrPopulation& mothers,
				const BitChrPopulation& fathers, std::uint_fast32_t seed0);
	std::size_t Morgan_to_index(double M) const;
	double get_length() const;
	std::vector<std::size_t> select_random_crossover_points(
												std::mt19937 &engine) const;
	void reduce(std::size_t parent_index,
				BitChrPopulation& new_population,
				std::size_t ind_index, std::size_t hap_index,
				std::mt19937 &engine) const;
	void write(std::ostream& os) const;
	BitChrPopulation *select(const std::vector<std::size_t>& indices) const;
	
public:
	static const BitChrPopulation *create_origins(std::size_t num_inds,
											const std::vector<GC::Pos>& ps,
											const ChromMap& cmap,
											const std::vector<double>& gratio,
											std::mt19937_64& engine);
	static std::vector<int> create_genotypes(std::size_t num_markers);
	static std::vector<GC::Pos> extract_positions_from_VCF(const VCF *vcf);
	static std::vector<uint64_t> create_genotypes_from_VCF(const VCF *vcf);
	static const BitChrPopulation *join(const BitChrPopulation *pop1,
										const BitChrPopulation *pop2);
};


//////////////////// Population ////////////////////

class Population {
public:
	using Pair = std::pair<std::size_t, std::size_t>;
	using Triplet = std::tuple<std::string, std::string, std::size_t>;
	
private:
	struct ConfigThread {
		const std::size_t	first;
		const std::size_t	num_threads;
		const Population&	mothers;
		const Population&	fathers;
		const std::vector<Pair>&	pairs;
		std::vector<BitChrPopulation *>&	chr_pops;
		const std::uint_fast32_t	seed0;
		
		ConfigThread(std::size_t i, std::size_t	nt,
						const Population& m, const Population& f,
						const std::vector<Pair>& ps,
						std::vector<BitChrPopulation *>& chr_ps,
						std::uint32_t s0) :
					first(i), num_threads(nt), mothers(m), fathers(f),
					pairs(ps), chr_pops(chr_ps), seed0(s0) { }
		
		std::size_t num_inds() const { return pairs.size(); }
		std::size_t num_chroms() const { return chr_pops.size(); }
	};
	
private:
	std::vector<const BitChrPopulation *>	chr_populations;
	const BaseInfo	*info;
	std::vector<std::string>	names;
	const std::vector<std::string>	mats;
	const std::vector<std::string>	pats;
	std::vector<std::vector<double>>	phenotypes;
	std::vector<const Trait *>	traits;
	
public:
	Population(const std::vector<const BitChrPopulation *>& chr_pops,
										const BaseInfo *bi,
										const std::vector<std::string>& ns,
										const std::vector<std::string>& ms,
										const std::vector<std::string>& ps) :
										chr_populations(chr_pops), info(bi),
										names(ns), mats(ms), pats(ps) { }
	~Population();
	
	const std::vector<std::string>& get_names() const { return names; }
	const std::vector<std::string>& get_mats() const { return mats; }
	const std::vector<std::string>& get_pats() const { return pats; }
	const std::string& get_name(std::size_t i) const { return names[i]; }
	const std::string& get_mat(std::size_t i) const { return mats[i]; }
	const std::string& get_pat(std::size_t i) const { return pats[i]; }
	std::size_t num_inds() const { return names.size(); }
	std::size_t num_chroms() const { return chr_populations.size(); }
	std::size_t num_markers() const;
	std::size_t num_traits() const;
	const BaseInfo *get_info() const { return info; }
	const ChromMap& get_chrmap(std::size_t i) const;
	const BitChrPopulation *get_chrpop(std::size_t i) const {
		return chr_populations[i];
	}
	const std::size_t get_chr_size(std::size_t i) const {
		return get_chrpop(i)->num_markers();
	}
	std::string get_genotype(std::size_t ind_index, std::size_t chr_index,
												std::size_t marker_id) const {
		return chr_populations[chr_index]->get_genotype(ind_index, marker_id);
	}
	// 0/0 => -1 0/1 => 0 1/1 => 1
	int get_int_genotype(std::size_t ind_index, std::size_t chr_index,
												std::size_t marker_id) const {
		return chr_populations[chr_index]->get_int_genotype(ind_index,
																marker_id);
	}
	int get_haplotype(std::size_t ind_index, std::size_t chr_index,
						std::size_t marker_index, std::size_t parent) const {
		return chr_populations[chr_index]->get_haplotype(ind_index,
														marker_index, parent);
	}
	void write_VCF(std::ostream& os) const;
	
	void set_phenotypes(const BaseInfo *info);
	void set_names(const std::vector<std::string>& names_);
	std::vector<std::vector<double>> compute_phenotypes(const BaseInfo *ifno,
														std::mt19937 &engine);
	std::vector<double> get_phenotypes(std::size_t i) const {
		return phenotypes[i];
	}
	void write_phenotypes(std::ostream& os) const;
	double mean(std::size_t i) const;
	double stddev(std::size_t i) const;
	void dispay_QTLs(std::size_t i) const;
	
	Population *select(const std::vector<std::size_t>& indices) const;
	
private:
	std::vector<double> select_phenotypes(
								const std::vector<std::size_t>& indices,
								size_t i) const;
	
public:
	static Population *create_origins(const BaseInfo *info,
										const std::vector<double>& gratio,
										const std::vector<std::string>& names);
	static Population *create_from_HaploArray(
								const std::vector<std::vector<uint64_t>>& genos,
								const BaseInfo *info,
								const std::vector<std::string>& names);
	static std::vector<Pair> make_pairs_randomly(std::size_t num_inds,
												const Population& mothers,
												const Population& fathers,
												std::mt19937& engine);
	static std::vector<Pair> make_pairs_by_table(
										const std::vector<Triplet>& table,
										const Population& mothers,
										const Population& fathers);
	static Population *cross_randomly(const Population& mothers,
								const Population& fathers, const BaseInfo *info,
								const std::vector<std::string>& names, int T);
	static Population *cross_by_table(const std::vector<Triplet>& table,
						const Population& mothers, const Population& fathers,
						const BaseInfo *info,
						const std::vector<std::string>& names, int T);
	static Population *cross(const std::vector<Pair>& pairs,
							const Population& mothers,
							const Population& fathers, const BaseInfo *info,
							const std::vector<std::string>& names,
							std::mt19937& engine, int T);
	static void cross_in_thread(void *config);
	static Population *join(const Population *pop1, const Population *pop2);
	static Population *create_from_VCF(const VCF *vcf, int seed);
	static std::size_t transform_NumericVector_to_geno(
									const Rcpp::NumericVector& haploArray,
									std::vector<uint64_t>& geno,
									std::size_t first, std::size_t last,
									std::size_t parent);
};

#endif

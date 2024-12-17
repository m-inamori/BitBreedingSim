#ifndef __POPULATION
#define __POPULATION

#include <ostream>
#include <vector>
#include <random>
#include "Map.h"
#include "Int.h"

class BaseInfo;


//////////////////// BitChrPopulation ////////////////////

class BitChrPopulation {
	using ConstIter = std::vector<Int::ull>::const_iterator;
	using Iter = std::vector<Int::ull>::iterator;
	using Pair = std::pair<std::size_t, std::size_t>;
	
private:
	std::vector<Int::ull>	genos;
	const std::size_t	num_inds;
	const ChromMap&	chrmap;
	
public:
	BitChrPopulation(const std::vector<Int::ull>& gs,
						std::size_t num, const ChromMap& cmap) :
									genos(gs), num_inds(num), chrmap(cmap) { }
	
	BitChrPopulation(std::size_t num_inds_, const ChromMap& cmap) :
							genos((cmap.get_num_markers()+63)/64*2*num_inds_),
							num_inds(num_inds_), chrmap(cmap) { }
	
	std::size_t num_markers() const { return chrmap.get_num_markers(); }
	std::size_t num_elements() const { return (num_markers()+63)/64; }
	std::size_t get_num_inds() const { return num_inds; }
	std::string get_genotype(std::size_t id_ind, std::size_t id_marker) const;
	int get_int_genotype(std::size_t id_ind, std::size_t id_marker) const;
	ConstIter get_haplotype(std::size_t ind_index, std::size_t hap_id) const {
		return genos.begin() + num_elements() * (ind_index * 2 + hap_id);
	}
	Iter get_mut_haplotype(std::size_t ind_index, std::size_t hap_id) {
		return genos.begin() + num_elements() * (ind_index * 2 + hap_id);
	}
	
	void cross(const std::vector<Pair>& pairs,
				const BitChrPopulation& mothers,
				const BitChrPopulation& fathers, std::uint_fast32_t seed0);
	void reduce(std::size_t parent_index,
				BitChrPopulation& new_population,
				std::size_t ind_index, std::size_t hap_index,
				std::mt19937 &engine) const;
	void write(std::ostream& os) const;
	BitChrPopulation *select(const std::vector<std::size_t>& indices) const;
	
public:
	static const BitChrPopulation *create_origins(std::size_t num_inds,
													const ChromMap& cmap,
													std::mt19937_64& engine);
	static std::vector<int> create_genotypes(std::size_t num_markers);
};


//////////////////// Population ////////////////////

class Population {
	using Pair = std::pair<std::size_t, std::size_t>;
	
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
	const Map& gmap;
	const std::vector<std::string>	names;
	std::vector<std::vector<double>>	phenotypes;
	std::vector<const Trait *>	traits;
	
public:
	Population(const std::vector<const BitChrPopulation *>& chr_pops,
					const Map& m, const std::vector<std::string>& ns) :
							chr_populations(chr_pops), gmap(m), names(ns) { }
	~Population();
	
	std::size_t num_inds() const { return names.size(); }
	std::size_t num_chroms() const { return chr_populations.size(); }
	const ChromMap&	get_chrmap(std::size_t i) const { return gmap.get_chr(i); }
	const BitChrPopulation	*get_chrpops(std::size_t i) const {
		return chr_populations[i];
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
	void write(std::ostream& os) const;
	
	void set_phenotypes(const BaseInfo *info);
	std::vector<std::vector<double>> compute_phenotypes(const BaseInfo *ifno,
														std::mt19937 &engine);
	std::vector<double> get_phenotypes(std::size_t i) const {
		return phenotypes[i];
	}
	void write_phenotypes(std::ostream& os) const;
	double mean(std::size_t i) const;
	double stddev(std::size_t i) const;
	void display_QTLs(std::size_t i) const;
	
	Population *select(const std::vector<std::size_t>& indices) const;
	
private:
	std::vector<double> select_phenotypes(
								const std::vector<std::size_t>& indices,
								size_t i) const;
	
public:
	static const Population *create_origins(std::size_t num_inds,
											const BaseInfo *info,
											const std::string& name_base);
	static Population *cross(std::size_t num_inds, const Population& mothers,
						const Population& fathers, const BaseInfo *info,
						const std::string& name_base, int T);
	static void cross_in_thread(void *config);
	static std::vector<Pair> make_pairs(std::size_t num_inds,
										const Population& mothers,
										const Population& fathers,
										std::mt19937& engine);
	static const BitChrPopulation *join(const BitChrPopulation *pop1,
										const BitChrPopulation *pop2);
};

#endif

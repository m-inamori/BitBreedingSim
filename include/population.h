#ifndef __POPULATION
#define __POPULATION

#include <vector>
#include <random>
#include "Map.h"


//////////////////// ChrPopulation ////////////////////

class ChrPopulation {
	using ConstIter = std::vector<int>::const_iterator;
	using Iter = std::vector<int>::iterator;
	using Pair = std::pair<std::size_t, std::size_t>;
	
	struct ConfigThread {
		const std::size_t	first;
		const std::size_t	num_threads;
		const std::vector<Pair>& pairs;
		const ChrPopulation&	mothers;
		const ChrPopulation&	fathers;
		ChrPopulation&	new_population;
		
		ConfigThread(std::size_t i, std::size_t	nt,
						const ChrPopulation& m, const ChrPopulation& f,
						const std::vector<Pair>& ps, ChrPopulation& new_pop) :
					first(i), num_threads(nt), pairs(ps),
					mothers(m), fathers(f), new_population(new_pop) { }
	};
	
private:
	std::vector<int>	genos;
	const ChromMap&	chrmap;
	
public:
	ChrPopulation(const std::vector<int>& gs, const ChromMap& cmap) :
												genos(gs), chrmap(cmap) { }
	
	ChrPopulation(std::size_t num_inds, const ChromMap& cmap) :
					genos(num_inds * cmap.num_markers() * 2), chrmap(cmap) { }
	
	std::size_t num_markers() const { return chrmap.num_markers(); }
	std::size_t num_inds() const { return genos.size() / (num_markers() * 2); }
	ConstIter get_haplotype(std::size_t ind_index, std::size_t hap_id) const {
		return genos.begin() + num_markers() * (ind_index * 2 + hap_id);
	}
	Iter get_mut_haplotype(std::size_t ind_index, std::size_t hap_id) {
		return genos.begin() + num_markers() * (ind_index * 2 + hap_id);
	}
	
	void cross(const std::vector<Pair>& pairs,
			const ChrPopulation& mothers, const ChrPopulation& fathers, int T);
	void cross_each(const ChrPopulation& mother,
					const ChrPopulation& father,
					std::size_t mat_index, std::size_t pat_index,
					std::size_t ind_index, std::mt19937 &engine);
	void reduce(std::size_t parent_index,
				ChrPopulation& new_population,
				std::size_t ind_index, std::size_t hap_index,
				std::mt19937 &engine) const;
	
public:
	static const ChrPopulation *create_origins(std::size_t num_inds,
												const ChromMap& cmap);
	static std::vector<int> create_genotypes(std::size_t num_markers);
	static void cross_in_thread(void *config);
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
		std::vector<ChrPopulation *>	chr_pops;
		
		ConfigThread(std::size_t i, std::size_t	nt,
						const Population& m, const Population& f,
						const std::vector<Pair>& ps) :
					first(i), num_threads(nt), mothers(m), fathers(f),
					pairs(ps), chr_pops(m.num_chroms()) { }
		
		std::size_t num_inds() const { return pairs.size(); }
		std::size_t num_chroms() const { return chr_pops.size(); }
	};
	
private:
	std::vector<const ChrPopulation *>	chr_populations;
	const Map& gmap;
	const std::vector<std::string>	names;
	
public:
	Population(const std::vector<const ChrPopulation *>& chr_pops,
				const Map& m, const std::vector<std::string>& ns) :
							chr_populations(chr_pops), gmap(m), names(ns) { }
	~Population();
	
	std::size_t num_inds() const { return names.size(); }
	std::size_t num_chroms() const { return chr_populations.size(); }
	const ChromMap&	get_chrmap(std::size_t i) const { return *gmap.get_chr(i); }
	const ChrPopulation	*get_chrpops(std::size_t i) const {
		return chr_populations[i];
	}
	
public:
	static const Population *create_origins(std::size_t num_inds,
											const Map& gmap,
											const std::string& name_base);
	static Population *cross(std::size_t num_inds,
						const Population& mothers, const Population& fathers,
						const Map& gmap, const std::string& name_base, int T);
	static void cross_in_thread(void *config);
	static std::vector<Pair> make_pairs(std::size_t num_inds,
										const Population& mothers,
										const Population& fathers);
};

#endif

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
	
	void reduce(std::size_t parent_index,
				ChrPopulation& new_population,
				std::size_t ind_index, std::size_t hap_index,
				std::mt19937 &engine) const;
	
public:
	static const ChrPopulation *create_origins(std::size_t num_inds,
												const ChromMap& cmap);
	static void cross(const std::vector<Pair>& pairs,
						const ChrPopulation& mathers,
						const ChrPopulation& fathers,
						ChrPopulation& new_population,
						std::mt19937 &engine);
	static std::vector<int> create_genotypes(std::size_t num_markers);
};


//////////////////// Population ////////////////////

class Population {
	using Pair = std::pair<std::size_t, std::size_t>;
	
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
	const ChrPopulation	*get_chrpops(std::size_t i) const {
		return chr_populations[i];
	}
	
public:
	static const Population *create_origins(std::size_t num_inds,
											const Map& gmap,
											const std::string& name_base);
	static Population *cross(std::size_t num_inds,
						const Population& mothers, const Population& fathers,
						const Map& gmap, const std::string& name_base);
	static std::vector<Pair> make_pairs(std::size_t num_inds,
										const Population& mothers,
										const Population& fathers);
};

#endif

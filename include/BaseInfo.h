#pragma once

#include <vector>
#include <string>
#include <random>

class Map;
class ChromMap;
class Trait;
class Population;


//////////////////// BaseInfo ////////////////////

class BaseInfo {
	const Map *gmap;
	std::vector<const Trait *>	traits;
	mutable std::mt19937	engine;
	
public:
	BaseInfo(const Map *m, const std::vector<const Trait *>& ts,
												std::uint_fast32_t seed);
	~BaseInfo();
	
	const Map& get_map() const { return *gmap; }
	std::size_t num_chroms() const;
	const ChromMap& get_chrom_map(std::size_t i) const;
	std::size_t num_traits() const { return traits.size(); }
	const Trait *get_trait(std::size_t i) const { return traits[i]; }
	const std::string& get_trait_name(std::size_t i) const;
	std::mt19937& get_random_engine() const { return engine; }
	
	std::vector<double> compute_phenotypes(const Population& pop,
											std::size_t trait_index) const;
	
	void set_trait();
	void set_trait_AD_multi(std::size_t num_loci, double h2, double H2);
	
public:
	static BaseInfo *create_default(int seed);
};

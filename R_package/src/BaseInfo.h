#pragma once

#include <vector>
#include <string>
#include <random>
#include "trait.h"

class Map;
class ChromMap;
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
	
	///// add A Trait /////
	void add_A_a_randomly(const std::string& name,
										const std::vector<Trait::Locus>& loci,
										double mean, double sd, double h2);
	void add_A_l_randomly(const std::string& name,
										const std::vector<double>& a,
										double mean, double h2);
	void add_A_al_randomly(const std::string& name,
										std::size_t num_loci,
										double mean, double sd, double h2);
	void add_A(const std::string& name, double mean, double h2,
										const std::vector<double>& a,
										const std::vector<Trait::Locus>& loci);
	
	///// add AD Trait /////
	void add_AD_a_randomly(const std::string& name,
										double mean, double h2, double H2,
										const std::vector<double>& ds,
										const std::vector<Trait::Locus>& loci);
	void add_AD_d_randomly(const std::string& name,
										double mean, double h2, double H2,
										const std::vector<double>& as,
										const std::vector<Trait::Locus>& loci);
	void add_AD_l_randomly(const std::string& name, double mean, double h2,
										const std::vector<double>& as,
										const std::vector<double>& ds);
	void add_AD_ad_randomly(const std::string& name, double mean, double sd,
										double h2, double H2,
										const std::vector<Trait::Locus>& loci);
	void add_AD_al_randomly(const std::string& name, double mean,
										double h2, double H2,
										const std::vector<double>& ds);
	void add_AD_dl_randomly(const std::string& name,
										double mean, double h2, double H2,
										const std::vector<double>& as);
	void add_AD_adl_randomly(const std::string& name, std::size_t num_loci,
										double mean, double sd,
										double h2, double H2);
	void add_AD(const std::string& name, double mean, double h2,
										const std::vector<double>& as,
										const std::vector<double>& ds,
										const std::vector<Trait::Locus>& loci);
	
public:
	static BaseInfo *create_default(int num_chroms, int num_markers,
										double cM, int bp, int seed);
};

#pragma once

#include <vector>
#include <string>
#include <random>
#include "trait.h"
#include "GenomicsCommon.h"

class Map;
class ChromMap;
class Population;

namespace GC = GenomicsCommon;


//////////////////// BaseInfo ////////////////////

class BaseInfo {
	const std::vector<std::vector<GC::Pos>>	positions;
	const Map *gmap;
	std::vector<const Trait *>	traits;
	mutable std::mt19937	engine;
	
public:
	BaseInfo(const std::vector<std::vector<GC::Pos>>& ps,
				const Map *m, const std::vector<const Trait *>& ts,
												std::uint_fast32_t seed) :
						positions(ps), gmap(m), traits(ts), engine(seed) { }
	~BaseInfo();
	
	const Map& get_map() const { return *gmap; }
	std::size_t num_chroms() const;
	const ChromMap& get_chrom_map(std::size_t i) const;
	std::size_t get_num_all_markers() const {
		std::size_t	num = 0;
		for(const auto& ps : positions) {
			num += ps.size();
		}
		return num;
	}
	std::size_t get_num_markers(std::size_t i) const {
		return positions[i].size();
	}
	const std::vector<GC::Pos>& get_positions(std::size_t i) const {
		return positions[i];
	}
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
	
	///// modify Trait /////
	void modify_trait_h2(std::size_t i, double h2);
	void modify_trait_h2_a(std::size_t i, double h2,
							const std::vector<double>& a);
	void modify_trait_h2_am(std::size_t i, double h2, double am);
	void modify_trait_a(std::size_t i, const std::vector<double>& a);
	void modify_trait_am(std::size_t i, double am);

public:
	static BaseInfo *create_default(
						const std::vector<std::vector<GC::Pos>>& positions,
						double cM, int seed);
	static BaseInfo *create_from_markers(
						const std::vector<std::size_t>& num_markers,
						const std::vector<int>& bps, std::uint_fast32_t seed);
};

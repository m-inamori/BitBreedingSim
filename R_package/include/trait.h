#pragma once

#include <vector>
#include <string>
#include <random>

#include "GenomicsCommon.h"
#include "exception_with_code.h"

class BaseInfo;
class Population;
class Map;

namespace GC = GenomicsCommon;


//////////////////// Trait ////////////////////

class Trait {
public:
	// (chrom index, marker index)
	using Locus = std::pair<std::size_t, std::size_t>;
	using Positions = std::vector<std::vector<GC::Pos>>;
	
private:
	const std::string	name;
	
public:
	Trait(const std::string& name_) : name(name_) { }
	virtual ~Trait() { }
	
	const std::string& get_name() const { return name; }
	virtual const std::string get_type() const = 0;
	virtual double phenotype(std::size_t ind_index, const Population& pop,
											std::mt19937& engine) const = 0;
	virtual double get_mean() const = 0;
	virtual double get_sd() const = 0;
	virtual double h2() const = 0;
	virtual double H2() const = 0;
	virtual std::size_t num_QTLs() const = 0;
	virtual std::vector<Locus> get_loci() const = 0;
	virtual std::vector<double> get_addivtives() const = 0;
	virtual std::vector<double> get_dominants() const = 0;
	virtual bool has_dominants() const = 0;
	
public:
	///// helper /////
	static Locus get_locus(std::size_t k, const Positions& positions);
	static std::size_t count_all_markers(const Positions& positions);
	static std::vector<double> decide_additives_randomly(
											std::size_t num_loci, double sd,
											double h2, std::mt19937& engine);
	static std::vector<double> decide_dominants_randomly(std::size_t num_loci,
												double sd, double h2, double H2,
												std::mt19937 &engine);
	static std::vector<Locus> decide_loci_randomly(std::size_t num_loci,
													const Positions& positions,
													std::mt19937 &engine);
	static double determine_sd_from_additives(const std::vector<double>& as,
														double h2);
	static double determine_sd_from_dominants(const std::vector<double>& ds,
														double h2, double H2);
	
	///// create A Trait /////
	static const Trait *create_A_a_randomly(const std::string& name,
										const std::vector<Trait::Locus>& loci,
										double mean, double sd, double h2,
										std::mt19937 &engine);
	static const Trait *create_A_l_randomly(const std::string& name,
											const std::vector<double>& a,
											double mean, double h2,
											const Positions& positions,
											std::mt19937 &engine);
	static const Trait *create_A_al_randomly(const std::string& name,
											 std::size_t num_loci,
											 double mean, double sd, double h2,
											 const Positions& positions,
											 std::mt19937 &engine);
	static const Trait *create_A(const std::string& name,
									double mean, double h2,
									const std::vector<double>& a,
									const std::vector<Locus>& loci);
	
	///// create AD Trait /////
	static const Trait *create_AD_a_randomly(const std::string& name,
											double mean, double h2, double H2,
											const std::vector<double>& ds,
											const std::vector<Locus>& loci,
											std::mt19937 &engine);
	static const Trait *create_AD_d_randomly(const std::string& name,
											double mean, double h2, double H2,
											const std::vector<double>& as,
											const std::vector<Locus>& loci,
											std::mt19937 &engine);
	static const Trait *create_AD_l_randomly(const std::string& name,
											double mean, double h2,
											const std::vector<double>& as,
											const std::vector<double>& ds,
											const Positions& positions,
											std::mt19937 &engine);
	static const Trait *create_AD_ad_randomly(const std::string& name,
											double mean, double sd,
											double h2, double H2,
											const std::vector<Locus>& loci,
											std::mt19937 &engine);
	static const Trait *create_AD_al_randomly(const std::string& name,
											double mean, double h2, double H2,
											const std::vector<double>& ds,
											const Positions& positions,
											std::mt19937 &engine);
	static const Trait *create_AD_dl_randomly(const std::string& name,
											double mean, double sd, double H2,
											const std::vector<double>& as,
											const Positions& positions,
											std::mt19937 &engine);
	static const Trait *create_AD_adl_randomly(const std::string& name,
											std::size_t num_loci,
											double mean, double sd,
											double h2, double H2,
											const Positions& positions,
											std::mt19937 &engine);
	static const Trait *create_AD(const std::string& name,
											double mean, double H2,
											const std::vector<double>& as,
											const std::vector<double>& ds,
											const std::vector<Locus>& loci);
};


//////////////////// TraitAOne ////////////////////

class TraitAOne : public Trait {
	const std::size_t	chr_index;
	const std::size_t	marker_index;
	const double	additive_effect;
	const double	mean;
	const double	error_std_dev;
	
public:
	TraitAOne(const std::string& name, std::size_t ci, std::size_t mi,
										double ae, double m, double esd) :
											Trait(name),
											chr_index(ci), marker_index(ci),
											additive_effect(ae), mean(m),
											error_std_dev(esd) { }
	~TraitAOne() { }
	
	const std::string get_type() const { return "Additive Effect Only"; }
	double phenotype(std::size_t ind_index, const Population& pop,
											std::mt19937& engine) const;
	double get_mean() const { return mean; }
	double get_sd() const;
	double h2() const;
	double H2() const { return h2(); }
	std::size_t num_QTLs() const { return 1; }
	std::vector<Locus> get_loci() const;
	std::vector<double> get_addivtives() const;
	std::vector<double> get_dominants() const;
	bool has_dominants() const { return false; }
};


//////////////////// TraitADOne ////////////////////

class TraitADOne : public Trait {
	const std::size_t	chr_index;
	const std::size_t	marker_index;
	const double	additive_effect;
	const double	dominant_effect;
	const double	mean;
	const double	error_std_dev;
	
public:
	TraitADOne(const std::string& name, std::size_t ci, std::size_t mi,
								double ae, double de, double m, double esd) :
											Trait(name),
											chr_index(ci), marker_index(ci),
											additive_effect(ae),
											dominant_effect(de), mean(m),
											error_std_dev(esd) { }
	~TraitADOne() { }
	
	const std::string get_type() const {
		return "Additive and Dominant Effect";
	}
	double phenotype(std::size_t ind_index, const Population& pop,
												std::mt19937& engine) const;
	double get_mean() const { return mean; }
	double get_sd() const;
	double h2() const;
	double H2() const;
	std::size_t num_QTLs() const { return 1; }
	std::vector<Locus> get_loci() const;
	std::vector<double> get_addivtives() const;
	std::vector<double> get_dominants() const;
	bool has_dominants() const { return true; }
};


//////////////////// TraitAMulti ////////////////////

class TraitAMulti : public Trait {
	const std::vector<Trait::Locus>	loci;
	const std::vector<double>	additive_effects;
	const double	mean;
	const double	error_std_dev;
	
public:
	TraitAMulti(const std::string& name,
				const std::vector<Trait::Locus>& ls,
				const std::vector<double>& ae, double m, double esd) :
					Trait(name), loci(ls),
					additive_effects(ae), mean(m), error_std_dev(esd) { }
	~TraitAMulti() { }
	
	const std::string get_type() const { return "Additive Effect Only"; }
	double phenotype(std::size_t ind_index, const Population& pop,
												std::mt19937& engine) const;
	double get_mean() const { return mean; }
	double get_sd() const;
	double h2() const;
	double H2() const { return h2(); }
	std::size_t num_loci() const { return loci.size(); }
	std::size_t num_QTLs() const { return loci.size(); }
	std::vector<Locus> get_loci() const;
	std::vector<double> get_addivtives() const;
	std::vector<double> get_dominants() const;
	bool has_dominants() const { return false; }
};


//////////////////// TraitADMulti ////////////////////

class TraitADMulti : public Trait {
	const std::vector<Trait::Locus>	loci;
	const std::vector<double>	additive_effects;
	const std::vector<double>	dominant_effects;
	const double	intercept;
	const double	error_std_dev;
	
public:
	TraitADMulti(const std::string& name,
					const std::vector<Trait::Locus>& ls,
					const std::vector<double>& ae,
					const std::vector<double>& de, double int_, double esd) :
								Trait(name), loci(ls),
								additive_effects(ae), dominant_effects(de),
								intercept(int_), error_std_dev(esd) { }
	~TraitADMulti() { }
	
	const std::string get_type() const {
		return "Additive and Dominant Effect";
	}
	double phenotype(std::size_t ind_index, const Population& pop,
												std::mt19937& engine) const;
	double get_mean() const;
	double get_sd() const;
	double h2() const;
	double H2() const;
	std::size_t num_loci() const { return loci.size(); }
	std::size_t num_QTLs() const { return loci.size(); }
	std::vector<Locus> get_loci() const;
	std::vector<double> get_addivtives() const;
	std::vector<double> get_dominants() const;
	bool has_dominants() const { return true; }
};

#pragma once

#include <vector>
#include <string>
#include <random>

#include "trait.h"

class Population;


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
	double calc_var() const;
	double get_sd() const;
	double h2() const;
	double H2() const { return h2(); }
	std::size_t num_QTLs() const { return 1; }
	std::vector<Locus> get_loci() const;
	std::vector<double> get_addivtives() const;
	std::vector<double> get_dominants() const;
	bool has_dominants() const { return false; }
	
	const Trait	*modify_h2_a(double h2, const std::vector<double>& a) const;
	const Trait *modify_h2_a_d(double h2, const std::vector<double>& a,
											const std::vector<double>& d) const;
	const Trait	*modify_a(const std::vector<double>& a) const;
	const Trait *modify_a_d(const std::vector<double>& a,
								const std::vector<double>& d) const;
};

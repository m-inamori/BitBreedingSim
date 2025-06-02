#pragma once

#include <vector>
#include <string>
#include <random>

#include "trait.h"

class Population;


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
	
	const Trait *modify(double h2_, double a, double d) const;
	const Trait	*modify_h2(double h2) const;
	const Trait	*modify_h2_a(double h2, const std::vector<double>& a) const;
	const Trait *modify_a(const std::vector<double>& a) const;
	const Trait *modify_h2_d(double h2, const std::vector<double>& d) const;
	const Trait *modify_h2_a_d(double h2, const std::vector<double>& a,
										const std::vector<double>& d) const;
	const Trait *modify_d(const std::vector<double>& d) const;
	const Trait *modify_a_d(const std::vector<double>& a,
										const std::vector<double>& d) const;
};

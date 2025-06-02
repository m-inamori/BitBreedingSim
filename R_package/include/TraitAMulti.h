#pragma once

#include <vector>
#include <string>
#include <random>

#include "trait.h"

class Population;


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
	
	const Trait	*modify(double h2, const std::vector<double>& a) const;
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

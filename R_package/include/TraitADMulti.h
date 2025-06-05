#pragma once

#include <vector>
#include <string>
#include <random>

#include "trait.h"

class Population;


//////////////////// TraitADMulti ////////////////////

class TraitADMulti : public Trait {
	const std::vector<Trait::Locus>	loci;
	const std::vector<double>	additive_effects;
	const std::vector<double>	dominant_effects;
	const double	mean;
	const double	error_std_dev;
	
public:
	TraitADMulti(const std::string& name,
					const std::vector<Trait::Locus>& ls,
					const std::vector<double>& ae,
					const std::vector<double>& de, double m, double esd) :
								Trait(name), loci(ls),
								additive_effects(ae), dominant_effects(de),
								mean(m), error_std_dev(esd) { }
	~TraitADMulti() { }
	
	const std::string get_type() const {
		return "Additive and Dominant Effect";
	}
	double phenotype(std::size_t ind_index, const Population& pop,
												std::mt19937& engine) const;
	double get_mean() const { return mean; }
	double calc_var() const;
	double get_sd() const;
	double h2() const;
	double H2() const;
	std::size_t num_loci() const { return loci.size(); }
	std::size_t num_QTLs() const { return loci.size(); }
	std::vector<Locus> get_loci() const { return loci; }
	std::vector<double> get_addivtives() const { return additive_effects; }
	std::vector<double> get_dominants() const { return dominant_effects; }
	bool has_dominants() const { return true; }
	
	const Trait *modify(const std::vector<double>& a,
						const std::vector<double>& d) const;
	const Trait	*modify_h2_a(double h2, const std::vector<double>& a) const;
	const Trait *modify_h2_a_d(double h2, const std::vector<double>& a,
											const std::vector<double>& d) const;
	const Trait *modify2(double h2_, const std::vector<double>& a,
										const std::vector<double>& d) const;
	const Trait	*modify_a(const std::vector<double>& a) const;
	const Trait *modify_a_d(const std::vector<double>& a,
								const std::vector<double>& d) const;
};

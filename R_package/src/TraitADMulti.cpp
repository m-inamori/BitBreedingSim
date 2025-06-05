#include <random>
#include <cmath>
#include "../include/TraitADMulti.h"
#include "../include/population.h"
#include "../include/common.h"

using namespace std;


//////////////////// TraitADMulti ////////////////////

double TraitADMulti::phenotype(size_t ind_index, const Population& pop,
												std::mt19937& engine) const {
	double	value = mean;
	for(size_t i = 0; i < num_loci(); ++i) {
		const size_t	chr_index = loci[i].first;
		const size_t	marker_index = loci[i].second;
		const int	gt = pop.get_int_genotype(ind_index, chr_index, marker_index);
		if(gt == 0)
			value += dominant_effects[i] / 2;
		else
			value += gt * additive_effects[i] - dominant_effects[i] / 2;
	}
	std::normal_distribution<> dist(0.0, error_std_dev);
	return value + dist(engine);
}

double TraitADMulti::calc_var() const {
	const double	ae2 = Common::dot_product(additive_effects,
												additive_effects);
	const double	de2 = Common::dot_product(dominant_effects,
												dominant_effects);
	return error_std_dev * error_std_dev + ae2 / 2 + de2 / 4;
}

double TraitADMulti::get_sd() const {
	return sqrt(calc_var());
}

double TraitADMulti::h2() const {
	const double	ae2 = Common::dot_product(additive_effects,
												additive_effects);
	const double	de2 = Common::dot_product(dominant_effects,
												dominant_effects);
	const double	sd2 = error_std_dev * error_std_dev + ae2 / 2 + de2 / 4;
	return ae2 / 2 / sd2;
}

double TraitADMulti::H2() const {
	const double	ae2 = Common::dot_product(additive_effects,
												additive_effects);
	const double	de2 = Common::dot_product(dominant_effects,
												dominant_effects);
	const double	sd2 = error_std_dev * error_std_dev + ae2 / 2 + de2 / 4;
	return (ae2 / 2 + de2 / 4) / sd2;
}

const Trait *TraitADMulti::modify(const vector<double>& a,
									const vector<double>& d) const {
	const double	a2 = Common::dot_product(a, a);
	const double	d2 = Common::dot_product(d, d);
	const double	new_esd = sqrt(calc_var() - a2 / 2 - d2 / 4);
	return new TraitADMulti(name, loci, a, d, mean, new_esd);
}

const Trait	*TraitADMulti::modify_a(const vector<double>& a) const {
	return modify(a, dominant_effects);
}

const Trait	*TraitADMulti::modify_a_d(const vector<double>& a,
									const vector<double>& d) const {
	return modify(a, d);
}

const Trait	*TraitADMulti::modify2(double h2_, const vector<double>& a,
											const vector<double>& d) const {
	const double	a2 = Common::dot_product(a, a);
	const double	d2 = Common::dot_product(d, d);
	const double	new_all_var = a2 / (h2_ * 2);
	const double	new_esd = sqrt(new_all_var - a2 / 2 - d2 / 4);
	return new TraitADMulti(name, loci, a, d, mean, new_esd);
}

const Trait	*TraitADMulti::modify_h2_a(double h2_,
										const vector<double>& a) const {
	return modify2(h2_, a, dominant_effects);
}

const Trait *TraitADMulti::modify_h2_a_d(double h2_,
											const vector<double>& a,
											const vector<double>& d) const {
	return modify2(h2_, a, d);
}

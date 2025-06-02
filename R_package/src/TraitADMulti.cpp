#include <random>
#include <cmath>
#include "../include/TraitADMulti.h"
#include "../include/population.h"
#include "../include/common.h"

using namespace std;


//////////////////// TraitADMulti ////////////////////

double TraitADMulti::phenotype(size_t ind_index, const Population& pop,
												std::mt19937& engine) const {
	double	value = intercept;
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

double TraitADMulti::get_mean() const {
	const double	sum_d = std::accumulate(dominant_effects.begin(),
											dominant_effects.end(), 0.0);
	return intercept + sum_d / 2;
}

double TraitADMulti::get_sd() const {
	const double	ae2 = Common::dot_product(additive_effects,
												additive_effects);
	const double	de2 = Common::dot_product(dominant_effects,
												dominant_effects);
	return sqrt(error_std_dev * error_std_dev + ae2 / 2 + de2 / 4);
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

const Trait *TraitADMulti::modify(double h2_, const vector<double>& a,
												const vector<double>& d) const {
	const double	ae2 = Common::dot_product(a, a);
	const double	de2 = Common::dot_product(d, d);
	const double	new_esd = sqrt((1.0-h2_)*ae2 - de2*h2_/2) / sqrt(2*h2_);
	return new TraitADMulti(name, loci, a, d, intercept, new_esd);
}

const Trait	*TraitADMulti::modify_h2(double h2_) const {
	return modify(h2_, additive_effects, dominant_effects);
}

const Trait	*TraitADMulti::modify_h2_a(double h2_,
										const vector<double>& a) const {
	return modify(h2_, a, dominant_effects);
}

const Trait *TraitADMulti::modify_a(const vector<double>& a) const {
	return modify(h2(), a, dominant_effects);
}

const Trait *TraitADMulti::modify_h2_d(double h2_,
										const vector<double>& d) const {
	return modify(h2_, additive_effects, d);
}

const Trait *TraitADMulti::modify_h2_a_d(double h2_,
											const vector<double>& a,
											const vector<double>& d) const {
	return modify(h2_, a, d);
}

const Trait *TraitADMulti::modify_d(const vector<double>& d) const {
	return modify(h2(), additive_effects, d);
}

const Trait *TraitADMulti::modify_a_d(const vector<double>& a,
										const vector<double>& d) const {
	return modify(h2(), a, d);
}

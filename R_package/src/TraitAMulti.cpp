#include <random>
#include <cmath>
#include "../include/TraitAMulti.h"
#include "../include/population.h"
#include "../include/common.h"

using namespace std;


//////////////////// TraitAMulti ////////////////////

double TraitAMulti::phenotype(size_t ind_index, const Population& pop,
												std::mt19937& engine) const {
	double	value = mean;
	for(size_t i = 0; i < num_loci(); ++i) {
		const size_t	chr_index = loci[i].first;
		const size_t	marker_index = loci[i].second;
		const int	gt = pop.get_int_genotype(ind_index,
												chr_index, marker_index);
		value += gt * additive_effects[i];
	}
	std::normal_distribution<> dist(0.0, error_std_dev);
	return value + dist(engine);
}

double TraitAMulti::get_sd() const {
	const double	ae2 = Common::dot_product(additive_effects,
												additive_effects);
	return sqrt(error_std_dev * error_std_dev + ae2 / 2);
}

double TraitAMulti::h2() const {
	const double	ae2 = Common::dot_product(additive_effects,
												additive_effects);
	const double	sd2 = error_std_dev * error_std_dev + ae2 / 2;
	return ae2 / 2 / sd2;
}

vector<Trait::Locus> TraitAMulti::get_loci() const {
	return loci;
}

vector<double> TraitAMulti::get_addivtives() const {
	return additive_effects;
}

vector<double> TraitAMulti::get_dominants() const {
	return vector<double>(num_QTLs(), 0.0);
}

const Trait	*TraitAMulti::modify(double h2_, const vector<double>& a) const {
	const double	ae2 = Common::dot_product(a, a);
	const double	new_esd = sqrt((1.0 - h2_) * ae2) / sqrt(2 * h2_);
	return new TraitAMulti(name, loci, a, mean, new_esd);
}

const Trait	*TraitAMulti::modify_h2(double h2_) const {
	return modify(h2_, additive_effects);
}

const Trait	*TraitAMulti::modify_h2_a(double h2_,
										const vector<double>& a) const {
	return modify(h2_, a);
}

const Trait *TraitAMulti::modify_a(const vector<double>& a) const {
	return modify(h2(), a);
}

const Trait *TraitAMulti::modify_h2_d(double h2_,
										const vector<double>& d) const {
	return modify(h2_, additive_effects);
}

const Trait *TraitAMulti::modify_h2_a_d(double h2_,
										const vector<double>& a,
										const vector<double>& d) const {
	return modify(h2_, a);
}

const Trait *TraitAMulti::modify_d(const vector<double>& d) const {
	return modify(h2(), additive_effects);
}

const Trait *TraitAMulti::modify_a_d(const vector<double>& a,
										const vector<double>& d) const {
	return modify(h2(), a);
}

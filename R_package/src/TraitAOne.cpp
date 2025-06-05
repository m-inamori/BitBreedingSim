#include <random>
#include <cmath>
#include "../include/TraitAOne.h"
#include "../include/population.h"

using namespace std;


//////////////////// TraitAOne ////////////////////

double TraitAOne::phenotype(size_t ind_index, const Population& pop,
												std::mt19937& engine) const {
	const int	gt = pop.get_int_genotype(chr_index, marker_index, ind_index);
	std::normal_distribution<> dist(0.0, error_std_dev);
	return gt * additive_effect + mean + dist(engine);
}

double TraitAOne::calc_var() const {
	return error_std_dev * error_std_dev +
			additive_effect * additive_effect / 2;
}

double TraitAOne::get_sd() const {
	return sqrt(calc_var());
}

double TraitAOne::h2() const {
	const double	ae2 = additive_effect * additive_effect;
	const double	sd2 = error_std_dev * error_std_dev + ae2 / 2;
	return ae2 / 2 / sd2;
}

vector<Trait::Locus> TraitAOne::get_loci() const {
	return vector<Locus>(1, Locus(chr_index, marker_index));
}

vector<double> TraitAOne::get_addivtives() const {
	return vector<double>(1, additive_effect);
}

vector<double> TraitAOne::get_dominants() const {
	return vector<double>(1, 0.0);
}

const Trait	*TraitAOne::modify_a(const vector<double>& a) const {
	const double	new_esd = sqrt(calc_var() - a[0] * a[0] / 2);
	return new TraitAOne(name, chr_index, marker_index, a[0], mean, new_esd);
}

const Trait *TraitAOne::modify_a_d(const vector<double>& a,
									const vector<double>& d) const {
	return modify_a(a);
}

const Trait	*TraitAOne::modify_h2_a(double h2_,
									const vector<double>& a) const {
	const double	new_all_var = a[0] * a[0] / (h2_ * 2);
	const double	new_esd = sqrt(new_all_var * (1.0 - h2_));
	return new TraitAOne(name, chr_index, marker_index, a[0], mean, new_esd);
}

const Trait *TraitAOne::modify_h2_a_d(double h2_,
										const vector<double>& a,
										const vector<double>& d) const {
	return modify_h2_a(h2_, a);
}

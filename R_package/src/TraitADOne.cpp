#include <random>
#include <cmath>
#include "../include/TraitADOne.h"
#include "../include/population.h"

using namespace std;


//////////////////// TraitADOne ////////////////////

double TraitADOne::phenotype(size_t ind_index, const Population& pop,
												std::mt19937& engine) const {
	const int	gt = pop.get_int_genotype(chr_index, marker_index, ind_index);
	std::normal_distribution<> dist(0.0, error_std_dev);
	if(gt == 0)
		return dominant_effect / 2 + mean + dist(engine);
	else
		return gt * additive_effect + mean - dominant_effect / 2 + dist(engine);
}

double TraitADOne::get_sd() const {
	return sqrt(error_std_dev * error_std_dev +
				additive_effect * additive_effect / 2 +
				dominant_effect * dominant_effect / 4);
}

double TraitADOne::h2() const {
	const double	ae2 = additive_effect * additive_effect;
	const double	de2 = dominant_effect * dominant_effect;
	const double	sd2 = error_std_dev * error_std_dev + ae2 / 2 + de2 / 4;
	return ae2 / 2 / sd2;
}

double TraitADOne::H2() const {
	const double	ae2 = additive_effect * additive_effect;
	const double	de2 = dominant_effect * dominant_effect;
	const double	sd2 = error_std_dev * error_std_dev + ae2 / 2 + de2 / 4;
	return (ae2 / 2 + de2 / 4) / sd2;
}

vector<Trait::Locus> TraitADOne::get_loci() const {
	return vector<Locus>(1, Locus(chr_index, marker_index));
}

vector<double> TraitADOne::get_addivtives() const {
	return vector<double>(1, additive_effect);
}

vector<double> TraitADOne::get_dominants() const {
	return vector<double>(1, dominant_effect);
}

const Trait *TraitADOne::modify(double h2_, double a, double d) const {
	const double	new_var = ((1.0-h2_)*a*a - d*d*h2_/2) / (2*h2_*h2_);
	const double	new_esd = sqrt(new_var);
	return new TraitADOne(name, chr_index, marker_index, a, d, mean, new_esd);
}

const Trait	*TraitADOne::modify_h2(double h2_) const {
	return modify(h2_, additive_effect, dominant_effect);
}

const Trait	*TraitADOne::modify_h2_a(double h2_,
										const vector<double>& a) const {
	return modify(h2_, a[0], dominant_effect);
}

const Trait *TraitADOne::modify_a(const vector<double>& a) const {
	return modify(h2(), a[0], dominant_effect);
}

const Trait *TraitADOne::modify_h2_d(double h2_,
										const vector<double>& d) const {
	return modify(h2_, additive_effect, d[0]);
}

const Trait *TraitADOne::modify_h2_a_d(double h2_,
										const vector<double>& a,
										const vector<double>& d) const {
	return modify(h2_, a[0], d[0]);
}

const Trait *TraitADOne::modify_d(const vector<double>& d) const {
	return modify(h2(), additive_effect, d[0]);
}

const Trait *TraitADOne::modify_a_d(const vector<double>& a,
									const vector<double>& d) const {
	return modify(h2(), a[0], d[0]);
}

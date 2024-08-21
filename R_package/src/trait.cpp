#include <random>
#include <cmath>
#include "trait.h"
#include "population.h"
#include "Map.h"

using namespace std;


//////////////////// Trait ////////////////////

const Trait *Trait::create_A(const string& name, double mean, double h2,
											const vector<double>& a,
											const vector<Locus>& loci) {
	// sd^2h2 = sum(a^2/4 + (-a)^2/4 for a in v) = sum(a^2/2 for a in v)
	// sd^2 = sum(a^2 for a in v)/2h2
	// error_var = s2^2(1-h2)
	double	s = 0.0;
	for(size_t i = 0; i < loci.size(); ++i)
		s += a[i] * a[i];
	const double	error_std_dev = sqrt(s * (1.0-h2) / (2.0-2*h2));
	if(a.size() == 1) {
		return new TraitAOne(name, loci[0].first, loci[0].second,
												a[0], mean, error_std_dev);
	}
	else {
		return new TraitAMulti(name, loci, a, mean, error_std_dev);
	}
}

const Trait *Trait::create_A_randomly(const string& name, size_t num_loci,
											double mean, double sd,
											double h2, const Map *gmap,
											std::mt19937 &engine) {
	const auto	a = decide_additives_randomly(sd, num_loci, h2, engine);
	const auto	loci = decide_loci_randomly(num_loci, gmap, engine);
	if(num_loci == 1) {
		return new TraitAOne(name, loci[0].first, loci[0].second,
											a[0], mean, sd * sqrt(1.0 - h2));
	}
	else {
		return new TraitAMulti(name, loci, a, mean, sd * sqrt(1.0 - h2));
	}
}

vector<double> Trait::decide_additives_randomly(size_t num_loci,
													double sd, double h2,
													std::mt19937 &engine) {
	if(num_loci == 1) {
		// if one locus, not random
		const double	a = sd * sqrt(2.0 * h2);
		return vector<double>(1, a);
	}
	else {
		// use random numbers to generate additive effects
		// generate additive effects according to a gamma distribution
		// and corrects for heritability
		// https://cpprefjp.github.io/reference/random/gamma_distribution.html
		std::gamma_distribution<>	geo_dist(1.0, 1.0);
		std::uniform_int_distribution<int>	dist_unif(0, 1);
		vector<double>	additive_effects(num_loci);
		double	s = 0.0;
		for(size_t i = 0; i < num_loci; ++i) {
			additive_effects[i] = geo_dist(engine) * (dist_unif(engine)*2-1);
			s += additive_effects[i] * additive_effects[i];
		}
		const double	c = sd * sqrt(2.0*h2/s);
		for(size_t i = 0; i < num_loci; ++i) {
			additive_effects[i] *= c;
		}
		return additive_effects;
	}
}

vector<Trait::Locus> Trait::decide_loci_randomly(size_t num_loci,
													const Map *gmap,
													std::mt19937 &engine) {
	vector<Locus>	loci(num_loci);
	const size_t	n = gmap->num_all_markers();
	std::uniform_int_distribution<>	dist(0, n - 1);
	for(size_t i = 0; i < num_loci; ++i) {
		const size_t	k = dist(engine);
		loci[i] = gmap->get_loci(k);
	}
	return loci;
}

const Trait *Trait::create_AD(const string& name, double mean, double h2,
												const vector<double>& as,
												const vector<double>& ds,
												const vector<Locus>& loci) {
	// sum_a2 = sum(a^2 for a in as)
	// sum_d2 = sum(d^2 for d in ds)
	// sum_d = sum(ds)
	// sd^2h2 = sum(a^2/4 + (-a)^2/4 for a in v) = sum_a2/2
	// sd^2 = sum_a2/2h2
	// sd^2H2 = sum_a2/2 + sum_d2/4
	// H2 = h2 + sum_d2/4sd^2
	// error_var = sd^2(1-H2)
	// mean = intercept + sum_d / 2
	double	sum_a2 = 0.0;
	double	sum_d2 = 0.0;
	double	sum_d = 0.0;
	for(size_t i = 0; i < loci.size(); ++i) {
		sum_a2 += as[i] * as[i];
		sum_d2 += ds[i] * ds[i];
		sum_d += ds[i];
	}
	const double	sd2 = sum_a2 / (2 * h2);
	const double	H2 = h2 + sum_d2 / (4 * sd2);
	const double	error_std_dev = sqrt(sd2 * (1.0 - H2));
	const double	intercept = mean - sum_d / 2;
	return new TraitADMulti(name, loci, as, ds, intercept, error_std_dev);
}

const Trait *Trait::create_AD_randomly(const string& name,
								double mean, double sd, double h2, double H2,
								const vector<Locus>& loci,
								std::mt19937 &engine) {
	const size_t	N = loci.size();
	if(N == 1) {
		// if one locus, not random
		const double	a = sd * sqrt(2.0 * h2);
		const double	error_std_dev = sd * sqrt(1.0 - h2);
		const double	d = sd * sqrt(4.0 * (H2 - h2));
		return new TraitADOne(name, loci[0].first, loci[0].second,
												a, d, mean, error_std_dev);
	}
	else {
		const auto	additive_effects = decide_additives_randomly(
														N, sd, h2, engine);
		
		// use normal distribution to generate dominant effects
		std::normal_distribution<>	norm_dist(0.0, sd * sqrt((H2 - h2) * 4 / N));
		vector<double>	dominant_effects(N);
		double	sum_d = 0.0;
		double	sum_d2 = 0.0;
		for(size_t i = 0; i < N; ++i) {
			dominant_effects[i] = norm_dist(engine);
			sum_d += dominant_effects[i];
			sum_d2 += dominant_effects[i] * dominant_effects[i];
		}
		// scaling for exact variance
		const double	c2 = sd * sqrt(4*(H2-h2)/sum_d2);
		for(size_t i = 0; i < N; ++i) {
			dominant_effects[i] *= c2;
		}
		
		const double	intercept = mean - c2 * sum_d / (2 * N);
		const double	error_std_dev = sd * sqrt(1.0 - H2);
		return new TraitADMulti(name, loci, additive_effects,
								dominant_effects, intercept, error_std_dev);
	}
}


//////////////////// TraitAOne ////////////////////

double TraitAOne::phenotype(size_t ind_index, const Population& pop,
												std::mt19937& engine) const {
	const int	gt = pop.get_int_genotype(chr_index, marker_index, ind_index);
	std::normal_distribution<> dist(0.0, error_std_dev);
	return gt * additive_effect + mean + dist(engine);
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

vector<Trait::Locus> TraitADOne::get_loci() const {
	return vector<Locus>(1, Locus(chr_index, marker_index));
}

vector<double> TraitADOne::get_addivtives() const {
	return vector<double>(1, additive_effect);
}

vector<double> TraitADOne::get_dominants() const {
	return vector<double>(1, dominant_effect);
}


//////////////////// TraitAMulti ////////////////////

double TraitAMulti::phenotype(size_t ind_index, const Population& pop,
												std::mt19937& engine) const {
	double	value = mean;
	for(size_t i = 0; i < num_loci(); ++i) {
		const size_t	chr_index = loci[i].first;
		const size_t	marker_index = loci[i].second;
		const int	gt = pop.get_int_genotype(chr_index, marker_index, ind_index);
		value += gt * additive_effects[i];
	}
	std::normal_distribution<> dist(0.0, error_std_dev);
	return value + dist(engine);
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

vector<Trait::Locus> TraitADMulti::get_loci() const {
	return loci;
}

vector<double> TraitADMulti::get_addivtives() const {
	return additive_effects;
}

vector<double> TraitADMulti::get_dominants() const {
	return dominant_effects;
}

#include <random>
#include <cmath>
#include "../include/trait.h"
#include "../include/BaseInfo.h"
#include "../include/population.h"
#include "../include/Map.h"
#include "../include/common.h"

using namespace std;


//////////////////// Trait ////////////////////

Trait::Locus Trait::get_locus(std::size_t k, const Positions& positions) {
	for(size_t i = 0; i < positions.size(); ++i) {
		if(k < positions[i].size()) {
			return Trait::Locus(i, k);
		}
		k -= positions[i].size();
	}
	return Trait::Locus(0, 0);	// not come here
}

size_t Trait::count_all_markers(const Positions& positions) {
	size_t	num = 0;
	for(const auto& p : positions) {
		num += p.size();
	}
	return num;
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

vector<double> Trait::decide_dominants_randomly(size_t num_loci,
												double sd, double h2, double H2,
												std::mt19937 &engine) {
	if(num_loci == 1) {
		// if one locus, not random
		const double	d = sd * sqrt(4.0 * (H2 - h2));
		return vector<double>(1, d);
	}
	else {
		// use normal distribution to generate dominant effects
		const double	sd_norm = sd * sqrt((H2 - h2) * 4 / num_loci);
		std::normal_distribution<>	norm_dist(0.0, sd_norm);
		vector<double>	ds(num_loci);
		for(size_t i = 0; i < num_loci; ++i) {
			ds[i] = norm_dist(engine);
		}
		const double	sum_d2 = Common::dot_product(ds, ds);
		// scaling for exact variance
		const double	c2 = sd * sqrt(4*(H2-h2)/sum_d2);
		for(size_t i = 0; i < num_loci; ++i) {
			ds[i] *= c2;
		}
		return ds;
	}
}

vector<Trait::Locus> Trait::decide_loci_randomly(size_t num_loci,
													const Positions& positions,
													std::mt19937 &engine) {
	vector<Locus>	loci(num_loci);
	
	const size_t	num = count_all_markers(positions);
	std::uniform_int_distribution<>	dist(0, num - 1);
	for(size_t i = 0; i < num_loci; ++i) {
		const size_t	k = dist(engine);
		loci[i] = Trait::get_locus(k, positions);
	}
	return loci;
}

double Trait::determine_sd_from_additives(const vector<double>& as, double h2) {
	const double	sum_a2 = Common::dot_product(as, as);
	const double	sd = sqrt(sum_a2 / (h2 * 2.0));
	return sd;
}

double Trait::determine_sd_from_dominants(const vector<double>& ds,
													double h2, double H2) {
	const double	sum_d2 = Common::dot_product(ds, ds);
	return sqrt(sum_d2 / ((H2 - h2) * 4.0));
}

const Trait *Trait::create_A_a_randomly(const string& name,
										const vector<Locus>& loci,
										double mean, double sd, double h2,
										std::mt19937 &engine) {
	const size_t	N = loci.size();
	const auto	as = decide_additives_randomly(N, sd, h2, engine);
	return create_A(name, mean, h2, as, loci);
}

const Trait *Trait::create_A_l_randomly(const string& name,
										const vector<double>& as,
										double mean, double h2,
										const Positions& positions,
										std::mt19937 &engine) {
	const size_t	N = as.size();
	const auto	loci = decide_loci_randomly(N, positions, engine);
	return create_A(name, mean, h2, as, loci);
}

const Trait *Trait::create_A_al_randomly(const string& name, std::size_t N,
										 double mean, double sd, double h2,
										 const Positions& positions,
										 std::mt19937 &engine) {
	const auto	as = decide_additives_randomly(N, sd, h2, engine);
	const auto	loci = decide_loci_randomly(N, positions, engine);
	return create_A(name, mean, h2, as, loci);
}

const Trait *Trait::create_A(const string& name, double mean, double h2,
											const vector<double>& as,
											const vector<Locus>& loci) {
	const double	sd = determine_sd_from_additives(as, h2);
	const double	error_std_dev = sd * sqrt(1.0 - h2);
	if(as.size() == 1) {
		return new TraitAOne(name, loci[0].first, loci[0].second,
												as[0], mean, error_std_dev);
	}
	else {
		return new TraitAMulti(name, loci, as, mean, error_std_dev);
	}
}

const Trait *Trait::create_AD_a_randomly(const string& name,
								double mean, double h2, double H2,
								const vector<double>& ds,
								const vector<Locus>& loci,
								std::mt19937 &engine) {
	const size_t	N = ds.size();
	const double	sd = determine_sd_from_dominants(ds, h2, H2);
	const auto	as = decide_additives_randomly(N, sd, h2, engine);
	return create_AD(name, mean, H2, as, ds, loci);
}

const Trait *Trait::create_AD_d_randomly(const string& name,
								double mean, double h2, double H2,
								const vector<double>& as,
								const vector<Locus>& loci,
								std::mt19937 &engine) {
	const size_t	N = as.size();
	const double	sd = determine_sd_from_additives(as, h2);
	const auto	ds = decide_dominants_randomly(N, sd, h2, H2, engine);
	return create_AD(name, mean, H2, as, ds, loci);
}

const Trait *Trait::create_AD_l_randomly(const string& name,
										 double mean, double h2,
										 const vector<double>& as,
										 const vector<double>& ds,
										 const Positions& positions,
										 std::mt19937 &engine) {
	const size_t	N = ds.size();
	const auto	loci = decide_loci_randomly(N, positions, engine);
	const double	sd = determine_sd_from_additives(as, h2);
	const double	sum_d2 = Common::dot_product(ds, ds);
	const double	H2 = h2 + sum_d2 / (4.0 * sd * sd);
	return create_AD(name, mean, H2, as, ds, loci);
}

const Trait *Trait::create_AD_ad_randomly(const string& name,
								double mean, double sd, double h2, double H2,
								const vector<Locus>& loci,
								std::mt19937 &engine) {
	const size_t	N = loci.size();
	const auto	as = decide_additives_randomly(N, sd, h2, engine);
	const auto	ds = decide_dominants_randomly(N, sd, h2, H2, engine);
	return create_AD(name, mean, H2, as, ds, loci);
}

const Trait *Trait::create_AD_al_randomly(const string& name,
										  double mean, double h2, double H2,
										  const vector<double>& ds,
										  const Positions& positions,
										  std::mt19937 &engine) {
	const size_t	N = ds.size();
	const double	sd = determine_sd_from_dominants(ds, h2, H2);
	const auto	as = decide_additives_randomly(N, sd, h2, engine);
	const auto	loci = decide_loci_randomly(N, positions, engine);
	return create_AD(name, mean, H2, as, ds, loci);
}

const Trait *Trait::create_AD_dl_randomly(const string& name,
										  double mean, double h2, double H2,
										  const vector<double>& as,
										  const Positions& positions,
										  std::mt19937 &engine) {
	const size_t	N = as.size();
	const double	sd = determine_sd_from_additives(as, h2);
	const vector<double> ds = decide_dominants_randomly(N, sd, h2, H2, engine);
	const auto	loci = decide_loci_randomly(N, positions, engine);
	return create_AD(name, mean, H2, as, ds, loci);
}

const Trait *Trait::create_AD_adl_randomly(const string& name, std::size_t N,
										   double mean, double sd,
										   double h2, double H2,
										   const Positions& positions,
										   std::mt19937 &engine) {
	const auto	as = decide_additives_randomly(N, sd, h2, engine);
	const auto	ds = decide_dominants_randomly(N, sd, h2, H2, engine);
	const auto	loci = decide_loci_randomly(N, positions, engine);
	return create_AD(name, mean, H2, as, ds, loci);
}

const Trait *Trait::create_AD(const string& name, double mean, double H2,
												const vector<double>& as,
												const vector<double>& ds,
												const vector<Locus>& loci) {
	const size_t	N = as.size();
	const double	sum_a2 = Common::dot_product(as, as);
	const double	sum_d2 = Common::dot_product(ds, ds);
	const double	sd = sqrt((sum_a2/2 + sum_d2/4) / H2);
	const double	error_std_dev = sd * sqrt(1.0 - H2);
	if(N == 1) {
		// if one locus, not random
		return new TraitADOne(name, loci[0].first, loci[0].second,
											as[0], ds[0], mean, error_std_dev);
	}
	else {
		double	sum_d = std::accumulate(ds.begin(), ds.end(), 0.0);
		const double	intercept = mean - sum_d / 2;
		return new TraitADMulti(name, loci, as, ds, intercept, error_std_dev);
	}
}


//////////////////// TraitAOne ////////////////////

double TraitAOne::phenotype(size_t ind_index, const Population& pop,
												std::mt19937& engine) const {
	const int	gt = pop.get_int_genotype(chr_index, marker_index, ind_index);
	std::normal_distribution<> dist(0.0, error_std_dev);
	return gt * additive_effect + mean + dist(engine);
}

double TraitAOne::get_sd() const {
	return sqrt(error_std_dev * error_std_dev +
				additive_effect * additive_effect / 2);
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

vector<Trait::Locus> TraitADMulti::get_loci() const {
	return loci;
}

vector<double> TraitADMulti::get_addivtives() const {
	return additive_effects;
}

vector<double> TraitADMulti::get_dominants() const {
	return dominant_effects;
}


//////////////////// Export ////////////////////


#include <random>
#include <cmath>
#include "../include/trait.h"
#include "../include/traitAOne.h"
#include "../include/traitADOne.h"
#include "../include/traitAMulti.h"
#include "../include/TraitADMulti.h"
#include "../include/BaseInfo.h"
#include "../include/population.h"
#include "../include/Map.h"
#include "../include/common.h"

using namespace std;
using namespace Rcpp;


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

// std::pairをRcppのリストに変換する関数
Rcpp::List wrap_pair(const Trait::Locus& p) {
	return Rcpp::List::create(_["first"] = p.first, _["second"] = p.second);
}

Rcpp::List Trait::get_info() const {
	// Convert Locus vector to Rcpp lists
	const std::vector<Trait::Locus>& loci = get_loci();
	Rcpp::List	loci_list(loci.size());
	for (std::size_t j = 0; j < loci.size(); ++j) {
		loci_list[j] = wrap_pair(loci[j]);
	}
	
	Rcpp::List	trait_list = Rcpp::List::create(
		_["name"] = get_name(),
		_["type"] = get_type(),
		_["mean"] = get_mean(),
		_["sd"] = get_sd(),
		_["h2"] = h2(),
		_["H2"] = H2(),
		_["loci"] = loci_list,
		_["additives"] = get_addivtives(),
		_["dominants"] = get_dominants(),
		_["hasdominants"] = has_dominants()
	);
	trait_list.attr("class") = "Trait";
	return trait_list;
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

#include <sstream>

#include "../include/BaseInfo.h"
#include "../include/Map.h"
#include "../include/population.h"
#include "../include/common.h"

using namespace std;
using namespace Rcpp;


//////////////////// BaseInfo ////////////////////

BaseInfo::~BaseInfo() {
	delete gmap;
	Common::delete_all(traits);
}

const ChromMap& BaseInfo::get_chrom_map(size_t i) const {
	return gmap->get_chr(i);
}

size_t BaseInfo::num_chroms() const {
	return gmap->num_chroms();
}

const string& BaseInfo::get_trait_name(size_t i) const {
	return traits[i]->get_name();
}

BaseInfo *BaseInfo::create_default(const vector<vector<GC::Pos>>& positions,
														double cM, int seed) {
	const size_t	num_chroms = positions.size();
	const Map	*gmap = Map::create_default(num_chroms, 1e8);
	vector<const Trait *>	traits;
	
	if(seed == -1) {
		std::random_device	seed_gen;
		return new BaseInfo(positions, gmap, traits, seed_gen());
	}
	else {
		const auto	s = static_cast<std::uint_fast32_t>(seed);
		return new BaseInfo(positions, gmap, traits, s);
	}
}

///// add A Trait /////
void BaseInfo::add_A_a_randomly(const string& name,
										const vector<Trait::Locus>& loci,
										double mean, double sd, double h2) {
	const Trait	*trait = Trait::create_A_a_randomly(name, loci,
													mean, sd, h2, engine);
	traits.push_back(trait);
}

void BaseInfo::add_A_l_randomly(const string& name,
										const vector<double>& a,
										double mean, double h2) {
	const Trait	*trait = Trait::create_A_l_randomly(name, a, mean, h2,
														positions, engine);
	traits.push_back(trait);
}

void BaseInfo::add_A_al_randomly(const string& name, size_t num_loci,
										double mean, double sd, double h2) {
	const Trait	*trait = Trait::create_A_al_randomly(name, num_loci, mean,
													sd, h2, positions, engine);
	traits.push_back(trait);
}

void BaseInfo::add_A(const string& name, double mean, double h2,
										const vector<double>& a,
										const vector<Trait::Locus>& loci) {
	const Trait	*trait = Trait::create_A(name, mean, h2, a, loci);
	traits.push_back(trait);
}

///// add AD Trait /////

void BaseInfo::add_AD_a_randomly(const string& name,
										double mean, double h2, double H2,
										const vector<double>& ds,
										const vector<Trait::Locus>& loci) {
	const Trait	*trait = Trait::create_AD_a_randomly(name, mean, h2, H2,
															ds, loci, engine);
	traits.push_back(trait);
}

void BaseInfo::add_AD_d_randomly(const string& name,
										double mean, double h2, double H2,
										const vector<double>& as,
										const vector<Trait::Locus>& loci) {
	const Trait	*trait = Trait::create_AD_d_randomly(name, mean, h2, H2,
															as, loci, engine);
	traits.push_back(trait);
}

void BaseInfo::add_AD_l_randomly(const string& name, double mean, double h2,
										const vector<double>& as,
										const vector<double>& ds) {
	const Trait	*trait = Trait::create_AD_l_randomly(name, mean, h2,
													as, ds, positions, engine);
	traits.push_back(trait);
}

void BaseInfo::add_AD_ad_randomly(const std::string& name, double mean,
										double sd, double h2, double H2,
										const vector<Trait::Locus>& loci) {
	const Trait	*trait = Trait::create_AD_ad_randomly(name, mean, sd, h2, H2,
																loci, engine);
	traits.push_back(trait);
}

void BaseInfo::add_AD_al_randomly(const std::string& name, double mean,
										double h2, double H2,
										const std::vector<double>& ds) {
	const Trait	*trait = Trait::create_AD_al_randomly(name, mean, h2, H2,
														ds, positions, engine);
	traits.push_back(trait);
}

void BaseInfo::add_AD_dl_randomly(const std::string& name,
										double mean, double h2, double H2,
										const std::vector<double>& as) {
	const Trait	*trait = Trait::create_AD_dl_randomly(name, mean, h2, H2,
														as, positions, engine);
	traits.push_back(trait);
}

void BaseInfo::add_AD_adl_randomly(const std::string& name,
										std::size_t num_loci, double mean,
										double sd, double h2, double H2) {
	const Trait	*trait = Trait::create_AD_adl_randomly(name, num_loci,
															mean, sd, h2, H2,
															positions, engine);
	traits.push_back(trait);
}

void BaseInfo::add_AD(const std::string& name, double mean, double H2,
										const std::vector<double>& as,
										const std::vector<double>& ds,
										const std::vector<Trait::Locus>& loci) {
	const Trait	*trait = Trait::create_AD(name, mean, H2, as, ds, loci);
	traits.push_back(trait);
}


///// misc /////

vector<double> BaseInfo::compute_phenotypes(const Population& pop,
												size_t trait_index) const {
	vector<double>	phenotypes(pop.num_inds());
	const Trait	*trait = traits[trait_index];
	for(size_t i = 0; i < pop.num_inds(); ++i) {
		phenotypes[i] = trait->phenotype(i, pop, engine);
	}
	return phenotypes;
}

// convert from an R data frame to a std::vector<Trait::Locus>
vector<Trait::Locus> dataframe_to_loci(const DataFrame& df) {
	vector<Trait::Locus> loci;
	IntegerVector chrom = df["chrom"];
	IntegerVector marker = df["marker"];
	for(int i = 0; i < chrom.size(); ++i) {
		loci.push_back(make_pair(chrom[i]-1, marker[i]-1));
	}
	return loci;
}

// [[Rcpp::export]]
SEXP createBaseInfoCpp(Rcpp::List pos, double cM, int seed) {
	vector<vector<GC::Pos>>	positions(pos.size());
	for(size_t i = 0; i < positions.size(); ++i) {
		NumericVector	r_vec = pos[i];
		positions[i].reserve(r_vec.size());
		for(double value : r_vec) {
			positions[i].push_back(static_cast<GC::Pos>(value));
		}
	}
	auto	*info = BaseInfo::create_default(positions, cM, seed);
	Rcpp::XPtr<BaseInfo> ptr(info, true);
	return ptr;
}

// [[Rcpp::export]]
SEXP createBaseInfoWithMap(Rcpp::List pos, Rcpp::List chrom_maps, int seed) {
	vector<vector<GC::Pos>>	positions(pos.size());
	for(size_t i = 0; i < positions.size(); ++i) {
		NumericVector	r_vec = pos[i];
		positions[i].reserve(r_vec.size());
		for(double value : r_vec) {
			positions[i].push_back(static_cast<GC::Pos>(value));
		}
	}
	const Map	*gmap = Map::create_map_from_list(chrom_maps);
	vector<const Trait*> traits;
	const auto	s = static_cast<std::uint_fast32_t>(seed);
	BaseInfo	*base_info = new BaseInfo(positions, gmap, traits, s);
	Rcpp::XPtr<BaseInfo> ptr(base_info, true);
	return ptr;
}

// [[Rcpp::export]]
int getNumChroms(SEXP ptr) {
	Rcpp::XPtr<BaseInfo> baseInfo(ptr);
	return baseInfo->num_chroms();
}

// [[Rcpp::export]]
int getNumTraits(SEXP ptr) {
	Rcpp::XPtr<BaseInfo> baseInfo(ptr);
	return baseInfo->num_traits();
}

// [[Rcpp::export]]
int getNumAllMarkers(SEXP ptr) {
	Rcpp::XPtr<BaseInfo> baseInfo(ptr);
	return baseInfo->get_num_all_markers();
}

// [[Rcpp::export]]
int getNumMarkers(SEXP ptr, std::size_t i) {
	Rcpp::XPtr<BaseInfo> baseInfo(ptr);
	return baseInfo->get_num_markers(i);
}

// [[Rcpp::export]]
SEXP getTraitCpp(SEXP baseInfoPtr, std::size_t i) {
	Rcpp::XPtr<BaseInfo>	ptr_info(baseInfoPtr);
	const Trait	*trait = ptr_info->get_trait(i);
	return trait->get_info();
}

// [[Rcpp::export]]
List getMapfromInfo(SEXP baseInfoPtr) {
	Rcpp::XPtr<BaseInfo>	ptr_info(baseInfoPtr);
	const Map&	gmap = ptr_info->get_map();
	
	// Convert Map to Rcpp List
	List	map_list(gmap.num_chroms());
	CharacterVector	names(gmap.num_chroms());
	for(size_t i = 0; i < gmap.num_chroms(); ++i) {
		const ChromMap&	cmap = gmap.get_chr(i);
		const auto	ps = cmap.collect_positions();
		const size_t	N = ps.size();
		NumericVector	cMs(N);
		IntegerVector	bps(N);
		for(size_t k = 0; k < N; ++k) {
			cMs[k] = ps[k].first;
			bps[k] = ps[k].second;
		}
		DataFrame	df = DataFrame::create(Named("cM") = cMs,
										   Named("position") = bps);
		
		map_list[i] = df;
		names[i] = cmap.get_name();
	}
	map_list.attr("names") = names;
	return map_list;
}

// [[Rcpp::export]]
void add_Trait_A_wrapper(SEXP ptr, std::string name, double mean, double h2,
						Nullable<double> sd_ = R_NilValue,
						Nullable<NumericVector> a = R_NilValue,
						Nullable<List> loci = R_NilValue, size_t num_loci = 1) {
	Rcpp::XPtr<BaseInfo> baseInfo(ptr);
	
	if(a.isNotNull() && loci.isNotNull()) {
		vector<double> a_vec = as<vector<double>>(a);
		vector<Trait::Locus> loci_vec = dataframe_to_loci(as<DataFrame>(loci));
		baseInfo->add_A(name, mean, h2, a_vec, loci_vec);
	}
	else if(a.isNotNull()) {
		vector<double> a_vec = as<vector<double>>(a);
		baseInfo->add_A_l_randomly(name, a_vec, mean, h2);
	}
	else if(loci.isNotNull()) {
		const double	sd = as<double>(sd_);
		vector<Trait::Locus> loci_vec = dataframe_to_loci(as<DataFrame>(loci));
		baseInfo->add_A_a_randomly(name, loci_vec, mean, sd, h2);
	}
	else {
		const double	sd = as<double>(sd_);
		baseInfo->add_A_al_randomly(name, num_loci, mean, sd, h2);
	}
}

// [[Rcpp::export]]
void add_Trait_AD_wrapper(SEXP baseInfoPtr, std::string name, double mean,
									Nullable<double> sd_ = R_NilValue,
									Nullable<double> h2_ = R_NilValue,
									Nullable<double> H2_ = R_NilValue,
									Nullable<NumericVector> a = R_NilValue,
									Nullable<NumericVector> ds = R_NilValue,
									Nullable<List> loci = R_NilValue,
									std::size_t num_loci = 1) {
	Rcpp::XPtr<BaseInfo> baseInfo(baseInfoPtr);
	
	if(loci.isNotNull()) {
		const double	H2 = as<double>(H2_);
		vector<Trait::Locus> loci_vec = dataframe_to_loci(as<DataFrame>(loci));
		if(a.isNotNull() && ds.isNotNull()) {
			std::vector<double> as_vec = as<std::vector<double>>(a);
			std::vector<double> ds_vec = as<std::vector<double>>(ds);
			baseInfo->add_AD(name, mean, H2, as_vec, ds_vec, loci_vec);
		}
		else if(a.isNotNull()) {
			const double	h2 = as<double>(h2_);
			std::vector<double> as_vec = as<std::vector<double>>(a);
			baseInfo->add_AD_d_randomly(name, mean, h2, H2, as_vec, loci_vec);
		}
		else if(ds.isNotNull()) {
			const double	h2 = as<double>(h2_);
			std::vector<double> ds_vec = as<std::vector<double>>(ds);
			baseInfo->add_AD_a_randomly(name, mean, h2, H2, ds_vec, loci_vec);
		}
		else {
			const double	sd = as<double>(sd_);
			const double	h2 = as<double>(h2_);
			baseInfo->add_AD_ad_randomly(name, mean, sd, h2, H2, loci_vec);
		}
	}
	else {
		if(a.isNotNull() && ds.isNotNull()) {
			const double	h2 = as<double>(h2_);
			std::vector<double> as_vec = as<std::vector<double>>(a);
			std::vector<double> ds_vec = as<std::vector<double>>(ds);
			baseInfo->add_AD_l_randomly(name, mean, h2, as_vec, ds_vec);
		}
		else if(ds.isNotNull()) {
			const double	h2 = as<double>(h2_);
			const double	H2 = as<double>(H2_);
			std::vector<double> ds_vec = as<std::vector<double>>(ds);
			baseInfo->add_AD_al_randomly(name, mean, h2, H2, ds_vec);
		}
		else if(a.isNotNull()) {
			const double	h2 = as<double>(h2_);
			const double	H2 = as<double>(H2_);
			std::vector<double> as_vec = as<std::vector<double>>(a);
			baseInfo->add_AD_dl_randomly(name, mean, h2, H2, as_vec);
		}
		else {
			const double	sd = as<double>(sd_);
			const double	h2 = as<double>(h2_);
			const double	H2 = as<double>(H2_);
			baseInfo->add_AD_adl_randomly(name, num_loci, mean, sd, h2, H2);
		}
	}
}

///// modify Trait /////

// [[Rcpp::export]]
void modify_Trait_Params_A_wrapper(SEXP baseInfoPtr, size_t i,
									Nullable<double> h2_ = R_NilValue,
									Nullable<NumericVector> a_ = R_NilValue,
									Nullable<double> am_ = R_NilValue) {
	Rcpp::XPtr<BaseInfo> info(baseInfoPtr);
	
	const Trait	*old_trait = info->get_trait(i-1);
	const auto	add = old_trait->get_addivtives();
	if(h2_.isNotNull()) {
		const double	h2 = as<double>(h2_);
		if(a_.isNull() && am_.isNull()) {
			info->set_trait(i-1, old_trait->modify_h2(h2));
		}
		else if(a_.isNotNull()) {
			const std::vector<double> a = as<std::vector<double>>(a_);
			info->set_trait(i-1, old_trait->modify_h2_a(h2, a));
		}
		else if(am_.isNotNull()) {
			const double	am = as<double>(am_);
			const auto	a = Common::multiply_by_constant(add, am);
			info->set_trait(i-1, old_trait->modify_h2_a(h2, a));
		}
	}
	else {
		if(a_.isNotNull()) {
			const std::vector<double> a = as<std::vector<double>>(a_);
			info->set_trait(i-1, old_trait->modify_a(a));
		}
		else if(am_.isNotNull()) {
			const double	am = as<double>(am_);
			const auto	a = Common::multiply_by_constant(add, am);
			info->set_trait(i-1, old_trait->modify_a(a));
		}
	}
	delete old_trait;
}

// [[Rcpp::export]]
void modify_Trait_Params_AD_wrapper(SEXP baseInfoPtr, size_t i,
									Nullable<double> h2_ = R_NilValue,
									Nullable<NumericVector> a_ = R_NilValue,
									Nullable<double> am_ = R_NilValue,
									Nullable<NumericVector> d_ = R_NilValue,
									Nullable<double> dm_ = R_NilValue) {
	Rcpp::XPtr<BaseInfo> info(baseInfoPtr);
	const Trait	*old_trait = info->get_trait(i-1);
	const auto	add = old_trait->get_addivtives();
	const auto	dom = old_trait->get_dominants();
	const Trait	*new_trait = nullptr;
	
	if(h2_.isNotNull()) {
		const double	h2 = as<double>(h2_);
		if(a_.isNull() && am_.isNull() && d_.isNull() && dm_.isNull()) {
			new_trait = old_trait->modify_h2(h2);
		}
		else if(a_.isNotNull() && d_.isNull() && dm_.isNull()) {
			const std::vector<double> a = as<std::vector<double>>(a_);
			new_trait = old_trait->modify_h2_a(h2, a);
		}
		else if(am_.isNotNull() && d_.isNull() && dm_.isNull()) {
			const double	am = as<double>(am_);
			const auto	a = Common::multiply_by_constant(add, am);
			new_trait = old_trait->modify_h2_a(h2, a);
		}
		else if(d_.isNotNull() && a_.isNull() && am_.isNull()) {
			const std::vector<double> d = as<std::vector<double>>(d_);
			new_trait = old_trait->modify_h2_d(h2, d);
		}
		else if(dm_.isNotNull() && a_.isNotNull() && dm_.isNotNull()) {
			const double	dm = as<double>(dm_);
			const auto	d = Common::multiply_by_constant(dom, dm);
			new_trait = old_trait->modify_h2_d(h2, d);
		}
		else if(a_.isNotNull() && d_.isNotNull()) {
			const std::vector<double> a = as<std::vector<double>>(a_);
			const std::vector<double> d = as<std::vector<double>>(d_);
			new_trait = old_trait->modify_h2_a_d(h2, a, d);
		}
		else if(a_.isNotNull() && dm_.isNotNull()) {
			const std::vector<double> a = as<std::vector<double>>(a_);
			const double	dm = as<double>(dm_);
			const auto	d = Common::multiply_by_constant(dom, dm);
			new_trait = old_trait->modify_h2_a_d(h2, a, d);
		}
		else if(am_.isNotNull() && d_.isNotNull()) {
			const double	am = as<double>(am_);
			const auto	a = Common::multiply_by_constant(add, am);
			const std::vector<double> d = as<std::vector<double>>(d_);
			new_trait = old_trait->modify_h2_a_d(h2, a, d);
		}
		else if(a_.isNotNull() && dm_.isNotNull()) {
			const double	am = as<double>(am_);
			const double	dm = as<double>(dm_);
			const auto	a = Common::multiply_by_constant(add, am);
			const auto	d = Common::multiply_by_constant(dom, dm);
			new_trait = old_trait->modify_h2_a_d(h2, a, d);
		}
	}
	else {
		if(a_.isNotNull() && d_.isNull() && dm_.isNull()) {
			const std::vector<double> a = as<std::vector<double>>(a_);
			new_trait = old_trait->modify_a(a);
		}
		else if(am_.isNotNull() && d_.isNull() && dm_.isNull()) {
			const double	am = as<double>(am_);
			const auto	a = Common::multiply_by_constant(add, am);
			new_trait = old_trait->modify_a(a);
		}
		else if(d_.isNotNull() && a_.isNull() && am_.isNull()) {
			const std::vector<double> d = as<std::vector<double>>(d_);
			new_trait = old_trait->modify_d(d);
		}
		else if(dm_.isNotNull() && a_.isNull() && am_.isNull()) {
			const double	dm = as<double>(dm_);
			const auto	d = Common::multiply_by_constant(dom, dm);
			new_trait = old_trait->modify_d(d);
		}
		else if(a_.isNotNull() && d_.isNotNull()) {
			const std::vector<double> a = as<std::vector<double>>(a_);
			const std::vector<double> d = as<std::vector<double>>(d_);
			new_trait = old_trait->modify_a_d(a, d);
		}
		else if(a_.isNotNull() && dm_.isNotNull()) {
			const std::vector<double> a = as<std::vector<double>>(a_);
			const double	dm = as<double>(dm_);
			const auto	d = Common::multiply_by_constant(dom, dm);
			new_trait = old_trait->modify_a_d(a, d);
		}
		else if(am_.isNotNull() && d_.isNotNull()) {
			const double	am = as<double>(am_);
			const auto	a = Common::multiply_by_constant(add, am);
			const std::vector<double> d = as<std::vector<double>>(d_);
			new_trait = old_trait->modify_a_d(a, d);
		}
		else if(am_.isNotNull() && dm_.isNotNull()) {
			const double	am = as<double>(am_);
			const double	dm = as<double>(dm_);
			const auto	a = Common::multiply_by_constant(add, am);
			const auto	d = Common::multiply_by_constant(dom, dm);
			new_trait = old_trait->modify_a_d(a, d);
		}
	}
	
	if(new_trait) {
		delete old_trait;
		info->set_trait(i-1, new_trait);
	}
}

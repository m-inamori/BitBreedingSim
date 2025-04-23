#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include "Map.h"
#include "common.h"

using namespace std;
using namespace Rcpp;


//////////////////// ChromMap ////////////////////

const ChromMap *ChromMap::create_default(const std::string& name,
														double slope) {
	return new ChromMapLinear(name, slope);
}


//////////////////// ChromMapLinear ////////////////////


//////////////////// ChromMapLines ////////////////////

GC::Pos ChromMapLines::Morgan_to_bp(double M) const {
	// Use binary search to find which line segment to use
	size_t	first = 0;
	size_t	last = Morgans.size();
	while(last - first > 1) {
		const size_t	mid = (first + last) / 2;
		if(Morgans[mid] < M)
			first = mid;
		else
			last = mid;
	}
	return (positions[first+1] - positions[first]) /
					(Morgans[first+1] - Morgans[first]) * (M - Morgans[first]) +
																Morgans[first];
}

double ChromMapLines::bp_to_Morgan(GC::Pos bp) const {
	size_t	first = 0;
	size_t	last = Morgans.size();
	while(last - first > 1) {
		const size_t	mid = (first + last) / 2;
		const double	m = positions[mid];
		if(Morgans[mid] < bp)
			last = mid;
		else
			first = mid;
	}
	return (Morgans[first+1] - Morgans[first]) /
			(positions[first+1] - positions[first]) * (bp - positions[first]) +
															positions[first];
}

ChromMapLines *ChromMapLines::create(const string& name,
									 const vector<pair<GC::Pos, double>>& gmap,
									 const vector<GC::Pos>& marker_positions) {
	// gmap[i]とgamp[i+1]で補完する
	auto	interpolate = [&gmap](size_t p, size_t i) {
		return gmap[i].second + (gmap[i+1].second - gmap[i].second) *
									(p - gmap[i].first) /
									(gmap[i+1].first - gmap[i].first);
	};
	
	vector<double>	Morgans;
	size_t	i = 0;	// index of gmap
	size_t	k = 0;	// index of marker_positions
	while(i < gmap.size() - 1 && k < marker_positions.size()) {
		const int	pos = marker_positions[k];
		if(pos < gmap[0].first) {
			Morgans.push_back(interpolate(pos, 0));
			++k;
		}
		else if(pos > gmap.back().first) {
			Morgans.push_back(interpolate(pos, gmap.size() - 2));
			++k;
		}
		else if(pos >= gmap[i+1].first) {
			++i;
		}
		else {
			Morgans.push_back(interpolate(pos, i));
			++k;
		}
	}
	
	return new ChromMapLines(name, Morgans, marker_positions);
}


//////////////////// Map ////////////////////

Map::~Map() {
	for(const auto& m: chr_maps) {
		delete m;
	}
}

size_t Map::num_points() const {
	size_t	num = 0;
	for(const auto& m : chr_maps) {
		num += m->num_points();
	}
	return num;
}

const Map *Map::create_default(size_t num_chroms, double slope) {
	vector<const ChromMap *>	maps(num_chroms);
	for(size_t i = 0; i < num_chroms; ++i) {
		stringstream	ss;
		ss << setfill('0') << setw(2) << "Chr" << i + 1;
		maps[i] = ChromMap::create_default(ss.str(), slope);
	}
	return new Map(maps);
}

bool Map::is_valid_map_line(const vector<string>& v) {
	return v.size() == 3 && Common::is_int(v[1]) && Common::is_double(v[2]);
}

vector<vector<string>> Map::read_lines(const string& path) {
	ifstream	ifs(path.c_str());
	if(!ifs)
		throw FileNotFoundException(path);
	
	vector<vector<string>>	table;
	vector<string>	not_three_columns_lines;
	string	line;
	while(getline(ifs, line)) {
		if(line.c_str()[line.length()-1] == '\r')
			line = line.substr(0, line.length()-1);
		istringstream	iss(line);
		vector<string>	v = Common::split(line, ',');
		if(is_valid_map_line(v) && not_three_columns_lines.size() < 5) {
			not_three_columns_lines.push_back(line);
		}
		table.push_back(v);
	}
	if(!not_three_columns_lines.empty()) {
		throw MapFormatException(not_three_columns_lines);
	}
	return table;
}

vector<pair<string, vector<pair<GC::Pos, double>>>>
			Map::divide_into_chromosomes(const vector<vector<string>>& table) {
	vector<pair<string, vector<pair<GC::Pos, double>>>>	gmaps;
	vector<pair<GC::Pos, double>>	gmap;
	string	chr;
	for(const auto& v : table) {
		if(v[0] != chr && !chr.empty()) {
			gmaps.push_back(make_pair(chr, gmap));
			gmap.clear();
			chr = v[0];
		}
		// guaranteed to convert to numbers
		const GC::Pos	pos = stoll(v[1]);
		const double	Morgan = stof(v[2]);
		gmap.push_back(make_pair(pos, Morgan));
	}
	gmaps.push_back(make_pair(chr, gmap));
	return gmaps;
}

const Map *Map::read(const string& path,
						const vector<GC::Pos>& marker_positions) {
	const auto	table = read_lines(path);
	const auto	maps = divide_into_chromosomes(table);
	vector<const ChromMap *>	chr_maps;
	for(const auto& m : maps) {
		chr_maps.push_back(ChromMapLines::create(m.first, m.second,
													marker_positions));
	}
	return new Map(chr_maps);
}

const ChromMap *Map::create_chrom_map_lines_from_df(Rcpp::DataFrame df,
													const string& name) {
	const vector<double>	cMs = Rcpp::as<vector<double>>(df["cM"]);
	const vector<GC::Pos>	positions = Rcpp::as<vector<GC::Pos>>(
															df["position"]);
	
	vector<double> morgans(cMs.size());
	std::transform(cMs.begin(), cMs.end(), morgans.begin(),
									[](double cm) { return cm / 100; });
	
	ChromMapLines* chrom_map = new ChromMapLines(name, morgans, positions);
	Rcpp::XPtr<ChromMapLines>	ptr(chrom_map, true);
	return ptr;
}

Map *Map::create_map_from_list(Rcpp::List chrom_maps) {
	std::vector<const ChromMap*> maps;
	Rcpp::CharacterVector names = chrom_maps.names();
	for(int i = 0; i < chrom_maps.size(); ++i) {
//		Rcpp::DataFrame	df = chrom_maps[i];
		Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(chrom_maps[i]);
		const string	name = Rcpp::as<string>(names[i]);
		const ChromMap	*chr_map = create_chrom_map_lines_from_df(df, name);
		maps.push_back(chr_map);
	}
	Map*	gmap = new Map(maps);
	return gmap;
}


//////////////////// MapFormatException ////////////////////

MapFormatException::MapFormatException(const vector<string>& lines) {
	stringstream	ss;
	if(lines.size() == 1)
		ss << "error : the following line doesn't have three columns :";
	else
		ss << "error : the following lines don't have three columns :";
	
	for(const auto& line : lines)
		ss << '\n' << line;
	
	message = ss.str();
}


//////////////////// Export ////////////////////

// [[Rcpp::export]]
SEXP getMapInfo(SEXP mapPtr) {
	Rcpp::XPtr<Map>	ptr_map(mapPtr);
	
	Rcpp::List	pop_list = Rcpp::List::create(
		_["num_chroms"] = ptr_map->num_chroms(),
		_["num_points"] = ptr_map->num_points()
	);
	return pop_list;
}

// [[Rcpp::export]]
SEXP getMapCpp(SEXP mapPtr) {
	Rcpp::XPtr<Map>	ptr_map(mapPtr);
	CharacterVector	chrs(ptr_map->num_chroms());
	IntegerVector	positions(ptr_map->num_chroms());
	size_t	k = 0;
	for(size_t i = 0; i < ptr_map->num_chroms(); ++i) {
		const ChromMap&	cmap = ptr_map->get_chr(i);
		const auto	ps = cmap.collect_positions();
		for(size_t j = 0; j < ps.size(); ++j) {
			chrs[k] = ptr_map->get_chr(i).get_name();
			positions[k] = ps[j].first;
			k += 1;
		}
	}
	
	Rcpp::List	pop_list = Rcpp::List::create(
		_["chrom"] = chrs,
		_["position"] = positions
	);
	return pop_list;
}

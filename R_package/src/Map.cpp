#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include "Map.h"
#include "common.h"

using namespace std;


//////////////////// ChromMap ////////////////////

size_t Map::num_all_markers() const {
	size_t	num = 0;
	for(const auto& m : chr_maps) {
		num += m->get_num_markers();
	}
	return num;
}

vector<size_t> ChromMap::select_random_crossover_points(
										std::mt19937 &engine) const {
	// where crossovers happen
	const double	L = get_length();
	vector<size_t>	pts;
	double	m0 = 0.0;	// start from 0.0 Morgan
	std::exponential_distribution<>	dist(1.0);
	while(true) {
		const double	m = m0 + dist(engine);
		if(m > L)
			break;
		pts.push_back(Morgan_to_index(m));
		m0 = m;
	}
	return pts;
}

const ChromMap *ChromMap::create_default(const std::string& name,
											std::size_t num_markers,
											double length, std::size_t bps) {
	return new ChromMapLinear(name, num_markers, length, bps);
}


//////////////////// ChromMapLinear ////////////////////


//////////////////// ChromMapLines ////////////////////

size_t ChromMapLines::Morgan_to_index(double M) const {
	size_t	first = 0;
	size_t	last = Morgans.size();
	while(first < last - 1) {
		const size_t	mid = (first + last) / 2;
		if(Morgans[mid] > M) {
			last = mid;
		}
		else {
			first = mid;
		}
	}
	
	if(M - Morgans[first] < Morgans[last] - M || last == Morgans.size()) {
		return first;
	}
	else {
		return last;
	}
}

ChromMapLines *ChromMapLines::create(const string& name,
									 const vector<pair<int, double>>& gmap,
									 const vector<int>& marker_positions) {
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

pair<size_t, size_t> Map::get_loci(size_t k) const {
	for(size_t i = 0; i < num_chroms(); ++i) {
		if(k < num_markers(i)) {
			return make_pair(i, k);
		}
		k -= num_markers(i);
	}
	return pair<size_t, size_t>(num_chroms(), 0);
}

const Map *Map::create_default(size_t num_chroms, size_t num_markers,
													double L, size_t bps) {
	vector<const ChromMap *>	maps(num_chroms);
	for(size_t i = 0; i < num_chroms; ++i) {
		stringstream	ss;
		ss << setfill('0') << setw(2) << "Chr" << i + 1;
		maps[i] = ChromMap::create_default(ss.str(), num_markers, L, bps);
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

vector<pair<string, vector<pair<int, double>>>> Map::divide_into_chromosomes(
										const vector<vector<string>>& table) {
	vector<pair<string, vector<pair<int, double>>>>	gmaps;
	vector<pair<int, double>>	gmap;
	string	chr;
	for(const auto& v : table) {
		if(v[0] != chr && !chr.empty()) {
			gmaps.push_back(make_pair(chr, gmap));
			gmap.clear();
			chr = v[0];
		}
		// guaranteed to convert to numbers
		const int	pos = stoi(v[1]);
		const double	Morgan = stof(v[2]);
		gmap.push_back(make_pair(pos, Morgan));
	}
	gmaps.push_back(make_pair(chr, gmap));
	return gmaps;
}

const Map *Map::read(const string& path, const vector<int>& marker_positions) {
	const auto	table = read_lines(path);
	const auto	maps = divide_into_chromosomes(table);
	vector<const ChromMap *>	chr_maps;
	for(const auto& m : maps) {
		chr_maps.push_back(ChromMapLines::create(m.first, m.second,
													marker_positions));
	}
	return new Map(chr_maps);
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

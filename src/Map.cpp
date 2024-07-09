#include <sstream>
#include <iomanip>
#include <random>
#include "../include/Map.h"

using namespace std;


//////////////////// ChromMap ////////////////////

vector<size_t> ChromMap::select_random_crossover_points(
										std::mt19937 &engine) const {
	// where crossovers happen
	const double	L = Morgans.back();
	vector<size_t>	pts;
	double	pt0 = 0.0;
	std::exponential_distribution<>	dist(1.0);
	while(true) {
		const double	pt = pt0 + dist(engine);
		if(pt > L)
			break;
		pts.push_back(static_cast<size_t>(pt * 1000.0));
		pt0 = pt;
	}
	return pts;
}

const ChromMap *ChromMap::create_default(const std::string& name,
												std::size_t N, double L) {
	vector<double>	ms(N);
	for(size_t i = 0; i < N; ++i) {
		ms[i] = (i + 1) * L / N;
	}
	return new ChromMap(name, ms);
}


//////////////////// Map ////////////////////

Map::~Map() {
	for(auto p = chr_maps.begin(); p != chr_maps.end(); ++p) {
		delete *p;
	}
}

const Map *Map::create_default(size_t num_chroms,
								size_t num_markers, double L) {
	vector<const ChromMap *>	maps(num_chroms);
	for(size_t i = 0; i < num_chroms; ++i) {
		stringstream	ss;
		ss << setfill('0') << setw(2) << "Chr" << i + 1;
		maps[i] = ChromMap::create_default(ss.str(), num_markers, L);
	}
	return new Map(maps);
}

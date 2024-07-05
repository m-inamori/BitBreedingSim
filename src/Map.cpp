#include <sstream>
#include <iomanip>
#include "../include/Map.h"

using namespace std;


//////////////////// ChromMap ////////////////////

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
		maps.push_back(ChromMap::create_default(ss.str(), num_markers, L));
	}
	return new Map(maps);
}

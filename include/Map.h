#ifndef __MAP
#define __MAP

#include <vector>
#include <string>
#include <random>


//////////////////// ChromMap ////////////////////

class ChromMap {
	const std::string	chr_name;
	const std::vector<double>	Morgans;	// each marker
	
public:
	ChromMap(const std::string& name,
			const std::vector<double>& ms) : chr_name(name), Morgans(ms) { }
	
	std::size_t num_markers() const { return Morgans.size(); }
	double length() const { return Morgans.back(); }
	std::vector<std::size_t> select_random_crossover_points(
								   std::mt19937 &engine) const;
	
public:
	static const ChromMap *create_default(const std::string& name,
										std::size_t num_markers, double length);
};


//////////////////// Map ////////////////////

class Map {
	const std::vector<const ChromMap *>	chr_maps;
	
public:
	Map(const std::vector<const ChromMap *>& ms) : chr_maps(ms) { }
	~Map();
	
	std::size_t num_chroms() const { return chr_maps.size(); }
	const ChromMap *get_chr(std::size_t i) const { return chr_maps[i]; }
	std::size_t num_markers(std::size_t i) const {
		return chr_maps[i]->num_markers();
	}
	
public:
	static const Map *create_default(std::size_t num_chroms,
										std::size_t num_markers, double L);
};

#endif

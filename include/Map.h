#ifndef __MAP
#define __MAP

#include <vector>
#include <string>


//////////////////// ChromMap ////////////////////

class ChromMap {
	const std::string	chr_name;
	const std::vector<double>	Morgans;	// each marker
	
public:
	ChromMap(const std::string& name,
			const std::vector<double>& ms) : chr_name(name), Morgans(ms) { }
	
	static const ChromMap *create_default(const std::string& name,
										std::size_t num_markers, double length);
};


//////////////////// Map ////////////////////

class Map {
	const std::vector<const ChromMap *>	chr_maps;
	
public:
	Map(const std::vector<const ChromMap *>& ms) : chr_maps(ms) { }
	~Map();
	
	const ChromMap *get_chr(std::size_t i) const { return chr_maps[i]; }
	
public:
	static const Map *create_default(std::size_t num_chroms,
										std::size_t num_markers, double L);
};

#endif

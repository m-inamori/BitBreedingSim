#pragma once

#include <vector>
#include <string>
#include <random>
#include <Rcpp.h>
#include "exception_with_code.h"
#include "GenomicsCommon.h"

namespace GC = GenomicsCommon;


//////////////////// ChromMap ////////////////////

class ChromMap {
	const std::string	chr_name;
	
public:
	ChromMap(const std::string& name) : chr_name(name) { }
	virtual ~ChromMap() { }
	
	virtual GC::Pos Morgan_to_bp(double M) const = 0;
	virtual double bp_to_Morgan(GC::Pos bp) const = 0;
	virtual std::vector<std::pair<GC::Pos, double>>
								collect_positions() const = 0;
	virtual std::size_t num_points() const = 0;
	
	const std::string& get_name() const { return chr_name; }
	
public:
	// slope: bp per Morgan
	static const ChromMap *create_default(const std::string& name,
														double slope);
};


//////////////////// ChromMapLinear ////////////////////

class ChromMapLinear : public ChromMap {
	const double	slope;	// bp per Morgan
	
public:
	ChromMapLinear(const std::string& name, double slope_) :
								ChromMap(name), slope(slope_) { }
	~ChromMapLinear() { }
	
	///// virtual method /////
	GC::Pos Morgan_to_bp(double M) const {
		return static_cast<int>(slope * M);
	}
	double bp_to_Morgan(GC::Pos bp) const { return slope * bp; }
	std::vector<std::pair<GC::Pos, double>> collect_positions() const {
		return std::vector<std::pair<GC::Pos, double>> {
			std::pair<GC::Pos, double>(0, 0.0),
			std::pair<GC::Pos, double>(static_cast<GC::Pos>(slope), 100.0)
		};
	}
	std::size_t num_points() const { return 2; }
};


//////////////////// ChromMapLines ////////////////////

// broken line
class ChromMapLines : public ChromMap {
	const std::vector<double>	Morgans;
	const std::vector<GC::Pos>	positions;
	
public:
	ChromMapLines(const std::string& name, const std::vector<double>& m,
											const std::vector<GC::Pos>& ps) :
								ChromMap(name), Morgans(m), positions(ps) { }
	~ChromMapLines() { }
	
	///// virtual method /////
	GC::Pos Morgan_to_bp(double M) const;
	double bp_to_Morgan(GC::Pos bp) const;
	std::vector<std::pair<GC::Pos, double>> collect_positions() const {
		std::vector<std::pair<GC::Pos, double>>	ps;
		for(std::size_t i = 0; i < positions.size(); ++i) {
			ps.push_back(std::pair<GC::Pos, double>(positions[i], Morgans[i]));
		}
		return ps;
	}
	std::size_t num_points() const { return Morgans.size(); }
	
public:
	static ChromMapLines *create(const std::string& name, 
							const std::vector<std::pair<GC::Pos, double>>& gmap,
							const std::vector<GC::Pos>& marker_positions);
};


//////////////////// Map ////////////////////

class Map {
	const std::vector<const ChromMap *>	chr_maps;
	
public:
	Map(const std::vector<const ChromMap *>& ms) : chr_maps(ms) { }
	~Map();
	
	std::size_t num_chroms() const { return chr_maps.size(); }
	const ChromMap& get_chr(std::size_t i) const { return *chr_maps[i]; }
	std::size_t num_points() const;
	
public:
	// slope: bp per Morgan
	static const Map *create_default(std::size_t num_chroms, double slope);
	static const Map *read(const std::string& path,
							const std::vector<GC::Pos>& marker_positions);
	static Map *create_map_from_list(Rcpp::List chrom_maps);

private:
	static std::vector<std::vector<std::string>>
						read_lines(const std::string& path);
	static bool is_valid_map_line(const std::vector<std::string>& v);
	static std::vector<std::pair<std::string,
								 std::vector<std::pair<GC::Pos, double>>>>
						divide_into_chromosomes(
							const std::vector<std::vector<std::string>>& table);
	static const ChromMap *create_chrom_map_lines_from_df(Rcpp::DataFrame df,
													const std::string& name);
};


//////////////////// MapFormatException ////////////////////

class MapFormatException : public ExceptionWithCode {
private:
    std::string	message;
	
public:
    MapFormatException(const std::vector<std::string>& lines);
    
    ErrorCode::Type get_error_code() const {
		return ErrorCode::MAP_INVALID_FORMAT;
	}
	
    const char *what() const noexcept override {
		return message.c_str();
	}
};

#pragma once

#include <vector>
#include <string>
#include <random>
#include <Rcpp.h>
#include "exception_with_code.h"


//////////////////// ChromMap ////////////////////

class ChromMap {
	const std::string	chr_name;
	
public:
	ChromMap(const std::string& name) : chr_name(name) { }
	virtual ~ChromMap() { }
	
	virtual std::size_t get_num_markers() const = 0;
	virtual double get_length() const = 0;
	virtual int get_position(std::size_t i) const = 0;
	virtual std::size_t Morgan_to_index(double M) const = 0;
	
	const std::string& get_name() const { return chr_name; }
	std::vector<std::size_t> select_random_crossover_points(
								   std::mt19937 &engine) const;
	
public:
	static const ChromMap *create_default(const std::string& name,
											std::size_t num_markers,
											double length, std::size_t bps);
};


//////////////////// ChromMapLinear ////////////////////

class ChromMapLinear : public ChromMap {
	const std::size_t	num_markers;
	const double	length;
	const std::size_t	bps;
	
public:
	ChromMapLinear(const std::string& name, std::size_t n,
										double l, std::size_t b) :
						ChromMap(name), num_markers(n), length(l), bps(b) { }
	~ChromMapLinear() { }
	
	///// virtual method /////
	std::size_t get_num_markers() const { return num_markers; }
	double get_length() const { return length; }
	int get_position(std::size_t i) const {
		return static_cast<int>(bps * (i + 1) / num_markers);
	}
	std::size_t Morgan_to_index(double M) const {
		return std::min(num_markers - 1,
						static_cast<std::size_t>(M / length * num_markers));
	}
};


//////////////////// ChromMapLines ////////////////////

class ChromMapLines : public ChromMap {
	const std::vector<double>	Morgans;
	const std::vector<int>		positions;
	
public:
	ChromMapLines(const std::string& name, const std::vector<double>& m,
											const std::vector<int>& ps) :
								ChromMap(name), Morgans(m), positions(ps) { }
	~ChromMapLines() { }
	
	///// virtual method /////
	std::size_t get_num_markers() const { return Morgans.size(); }
	double get_length() const { return Morgans.back(); }
	int get_position(std::size_t i) const { return positions[i]; }
	std::size_t Morgan_to_index(double M) const;
	
public:
	static ChromMapLines *create(const std::string& name, 
								const std::vector<std::pair<int, double>>& gmap,
								const std::vector<int>& marker_positions);
};


//////////////////// Map ////////////////////

class Map {
	const std::vector<const ChromMap *>	chr_maps;
	
public:
	Map(const std::vector<const ChromMap *>& ms) : chr_maps(ms) { }
	~Map();
	
	std::size_t num_chroms() const { return chr_maps.size(); }
	const ChromMap& get_chr(std::size_t i) const { return *chr_maps[i]; }
	std::size_t num_markers(std::size_t i) const {
		return chr_maps[i]->get_num_markers();
	}
	std::size_t num_all_markers() const;
	std::pair<std::size_t, std::size_t> get_loci(std::size_t k) const;
	
public:
	static const Map *create_default(std::size_t num_chroms,
										std::size_t num_markers,
										double L, std::size_t bps);
	static const Map *read(const std::string& path,
							const std::vector<int>& marker_positions);
	static Map *create_map_from_list(Rcpp::List chrom_maps);

private:
	static std::vector<std::vector<std::string>>
						read_lines(const std::string& path);
	static bool is_valid_map_line(const std::vector<std::string>& v);
	static std::vector<std::pair<std::string, std::vector<std::pair<int, double>>>>
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

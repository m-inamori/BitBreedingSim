#ifndef __VCF
#define __VCF

#include <iostream>
#include <vector>
#include <string>
#include "filereader.h"
#include "GenomicsCommon.h"

namespace GC = GenomicsCommon;

using STRVEC = std::vector<std::string>;


//////////////////// VCFReader ////////////////////

class VCFReader {
	FileReaderBase	*reader;
	std::vector<STRVEC>	header;
	
public:
	explicit VCFReader(const std::string& path);
	VCFReader(const VCFReader&) = delete;
	VCFReader& operator=(const VCFReader&) = delete;
	~VCFReader();
	
	void read_header();
	const std::vector<STRVEC>& get_header() { return header; }
	STRVEC get_samples() const;
	STRVEC next();
};


//////////////////// VCF ////////////////////

class VCF {
	std::vector<STRVEC>	header;
	STRVEC				samples;
	std::vector<STRVEC>	data;
	
public:
	VCF(const std::vector<STRVEC>& h, const STRVEC& s,
			const std::vector<STRVEC>& d) : header(h), samples(s), data(d) { }
	~VCF() { }
	
	const std::vector<STRVEC>& get_header() const { return header; }
	const STRVEC& get_samples() const { return samples; }
	std::size_t num_samples() const { return samples.size(); }
	std::size_t size() const { return data.size(); }
	const std::vector<STRVEC>& get_data() const { return data; }
	const GC::Pos pos(std::size_t) const;
	const std::string& get_gt(std::size_t i, std::size_t k) const {
		return data[i][k+9];
	}
	const std::string& chrom(std::size_t i) const { return data[i][0]; }
	
	static VCF *read(const std::string& path);
	
	static void write_header(std::ostream& os,
							const std::vector<std::string>& samples);
	static void write_data_line(std::ostream& os, const std::string& chr,
							int pos, const std::vector<std::string>& gts);
};


//////////////////// VCFDivisor ////////////////////

class VCFDivisor {
	const VCF&	vcf;
	std::size_t	pos;
	
public:
	explicit VCFDivisor(const VCF& vcf_) : vcf(vcf_), pos(0) { }
	VCFDivisor(const VCFDivisor&) = delete;
	VCFDivisor& operator=(const VCFDivisor&) = delete;
	~VCFDivisor() { }
	
	const std::vector<STRVEC>& get_header() { return vcf.get_header(); }
	VCF *next();
};
#endif

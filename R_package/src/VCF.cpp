#include "../include/VCF.h"
#include "../include/population.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFReader ////////////////////

VCFReader::VCFReader(const string& path) {
	if(path.substr(path.length() - 3) != ".gz")
		reader = new FileReader(path);
	else
		reader = new FileReaderGZ(path);
}

VCFReader::~VCFReader() {
	delete reader;
}

void VCFReader::read_header() {
	string	line;
	while(reader->getline(line)) {
		header.push_back(Common::split(line, '\t'));
		if(line.substr(0, 6) == "#CHROM")
			break;
	}
}

STRVEC VCFReader::next() {
	string	line;
	if(!reader->getline(line))
		return STRVEC();
	
	return Common::split(line, '\t');
}

STRVEC VCFReader::get_samples() const {
	const STRVEC&	v = header.back();
	STRVEC	samples(v.begin() + 9, v.end());
	return samples;
}


//////////////////// VCF ////////////////////

const GC::Pos VCF::pos(std::size_t i) const {
	return stoll(data[i][1]);
}

VCF *VCF::read(const string& path) {
	VCFReader	*reader = new VCFReader(path);
	reader->read_header();
	const vector<STRVEC>	header = reader->get_header();
	const STRVEC	samples(header.back().begin() + 9, header.back().end());
	
	vector<STRVEC>	data;
	while(true) {
		const STRVEC	v = reader->next();
		if(v.empty())
			break;
		data.push_back(v);
	}
	delete reader;
	
	return new VCF(header, samples, data);
}

// add pop and create new vcf object
VCF *VCF::add_pop(const Population *pop) const {
	STRVEC	new_samples = samples;
	const STRVEC&	names = pop->get_names();
	new_samples.insert(new_samples.end(), names.begin(), names.end());
	
	vector<STRVEC>	new_header = header;
	STRVEC&	line = new_header.back();
	line.insert(line.end(), names.begin(), names.end());
	
	vector<STRVEC>	new_data = data;
	size_t	i = 0;	// marker serial number
	for(size_t chr_id = 0; chr_id < pop->num_chroms(); ++chr_id) {
		const auto	*chr_pop = pop->get_chrpop(chr_id);
		for(size_t j = 0; j < chr_pop->num_markers(); ++j) {
			for(size_t k = 0; k < names.size(); ++k) {
				new_data[i].push_back(chr_pop->get_genotype(k, j));
			}
			++i;
		}
	}
	
	return new VCF(new_header, new_samples, new_data);
}


void VCF::write_default_header(ostream& os, const STRVEC& samples) {
	os << "##fileformat=VCFv4.2\n";
	os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	os << "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT";
	for(const auto& sample : samples)
		os << '\t' << sample;
	os << "\n";
}

void VCF::write_data_line(std::ostream& os, const std::string& chr,
								int pos, const std::vector<std::string>& gts) {
	os << chr << '\t' << pos << '\t' << '.'
			<< '\t' << 'A' << '\t' << 'C' << '\t' << '.'
			<< '\t' << "PASS" << '\t' << '.' << '\t' << "GT";
	for(const auto& gt : gts)
		os << '\t' << gt;
	os << "\n";
}

void VCF::write_header(ostream& os) const {
	for(auto p = header.begin(); p != header.end(); ++p) {
		Common::write_tsv(*p, os);
	}
}

void VCF::write(ostream& os) const {
	write_header(os);
	for(auto p = data.begin(); p != data.end(); ++p) {
		Common::write_tsv(*p, os);
	}
}


//////////////////// VCFDivisor ////////////////////

VCF *VCFDivisor::next() {
	if(pos >= vcf.size())
		return NULL;
	
	const string&	prev_chr = vcf.chrom(pos);
	const vector<STRVEC>&	data = vcf.get_data();
	for(size_t new_pos = pos + 1; new_pos < vcf.size(); ++new_pos) {
		if(vcf.chrom(new_pos) != prev_chr) {
			vector<STRVEC>	new_data(data.begin() + pos,
											data.begin() + new_pos);
			pos = new_pos;
			return new VCF(vcf.get_header(), vcf.get_samples(), new_data);
		}
	}
	{
		vector<STRVEC>	new_data(data.begin() + pos, data.end());
		pos = data.size();
		return new VCF(vcf.get_header(), vcf.get_samples(), new_data);
	}
}

// [[Rcpp::export]]
SEXP readVCF(const std::string& filename) {
	Rcpp::XPtr<VCF> ptr(const_cast<VCF *>(VCF::read(filename)));
	return ptr;
}


//////////////////// Export ////////////////////

// [[Rcpp::export]]
int getNumIndsVCF(SEXP vcf) {
	Rcpp::XPtr<VCF> vcf_cpp(vcf);
	return vcf_cpp.get()->num_samples();
}

// [[Rcpp::export]]
int getNumMarkersVCF(SEXP vcf) {
	Rcpp::XPtr<VCF> vcf_cpp(vcf);
	return vcf_cpp.get()->size();
}

// [[Rcpp::export]]
SEXP addPopToVCF(SEXP vcf, SEXP pop) {
	Rcpp::XPtr<VCF> vcf_cpp(vcf);
	Rcpp::XPtr<Population> pop_cpp(pop);
	auto	*new_vcf = vcf_cpp.get()->add_pop(pop_cpp);
	Rcpp::XPtr<VCF> ptr(new_vcf, true);
	ptr.attr("class") = "VCF";
	return ptr;
}

// [[Rcpp::export]]
void writeVCF(SEXP vcf, const std::string& filename) {
	Rcpp::XPtr<VCF> vcfCpp(vcf);
	std::ofstream	ofs(filename);
	if(!ofs.is_open()) {
		Rcpp::stop("Failed to open file: " + filename);
	}
	vcfCpp.get()->write(ofs);
	ofs.close();
}

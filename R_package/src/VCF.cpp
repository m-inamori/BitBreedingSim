#include "VCF.h"
#include "common.h"

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


void VCF::write_header(ostream& os, const STRVEC& samples) {
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

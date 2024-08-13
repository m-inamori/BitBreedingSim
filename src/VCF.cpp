#include "../include/VCF.h"

void VCF::write_header(std::ostream& os,
						const std::vector<std::string>& samples) {
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

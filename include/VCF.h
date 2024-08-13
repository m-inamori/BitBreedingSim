#ifndef __VCF
#define __VCF

#include <iostream>
#include <vector>
#include <string>


//////////////////// VCF ////////////////////

namespace VCF {
	void write_header(std::ostream& os,
							const std::vector<std::string>& samples);
	void write_data_line(std::ostream& os, const std::string& chr,
							int pos, const std::vector<std::string>& gts);
}

#endif

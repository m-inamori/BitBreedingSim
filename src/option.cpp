#include <iostream>
#include <stdexcept>
#include "../include/option.h"

using namespace std;


//////////////////// Option ////////////////////

const Option *Option::create(int argc, char **argv) {
	if(argc != 3 && argc != 5)
		return nullptr;
	
	try {
		const size_t	ni = std::stoull(argv[1]);
		const size_t	nc = std::stoull(argv[2]);
		const string	str_num = flag_value("-t", argc, argv);
		const int		T = str_num.empty() ? 1 : std::stoi(str_num);
		return new Option(ni, nc, T);
	}
	catch(const invalid_argument& e) {
		return nullptr;
	}
}

string Option::flag_value(const string& s, int argc, char **argv) {
	for(size_t i = 1; i < (size_t)(argc - 1); ++i) {
		if(argv[i] == s)
			return argv[i+1];
	}
	return string();
}

void Option::usage() {
	cerr << "usage : BreedingSim num_inds num_chroms." << endl;
}

#include <iostream>
#include <stdexcept>
#include "../include/option.h"

using namespace std;


//////////////////// Option ////////////////////

const Option *Option::create(int argc, char **argv) {
	if(argc != 3 && argc != 5 && argc != 7)
		return nullptr;
	
	try {
		const size_t	ni = std::stoull(argv[1]);
		const size_t	nc = std::stoull(argv[2]);
		const int		T = flag_value_int("-t", argc, argv, 1);
		const string	path_out = flag_value("-o", argc, argv);
		return new Option(ni, nc, T, path_out);
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

int Option::flag_value_int(const string& s, int argc, char **argv,
													int default_value) {
	const string	str_value = flag_value(s, argc, argv);
	if(str_value.empty())
		return default_value;
	else
		return std::stoi(str_value);
}

void Option::usage() {
	cerr << "usage : BreedingSim num_inds num_chroms [-t num threads] [-o path out]." << endl;
}

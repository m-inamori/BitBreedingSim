#include <iostream>
#include <stdexcept>
#include "../include/option.h"

using namespace std;


//////////////////// Option ////////////////////

const Option *Option::create(int argc, char **argv) {
	if(!(argc % 2 == 1 && 3 <= argc && argc <= 11))
		return nullptr;
	
	try {
		const size_t	ni = std::stoull(argv[1]);
		const size_t	nc = std::stoull(argv[2]);
		const int		seed = flag_value_seed("-s", argc, argv);
		const int		T = flag_value_int("-t", argc, argv, 1);
		const string	path_out = flag_value("-o", argc, argv);
		const string	pheno_out = flag_value("-p", argc, argv);
		return new Option(ni, nc, seed, T, path_out, pheno_out);
	}
	catch(const invalid_argument& e) {
		cerr << e.what() << endl;
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

int Option::flag_value_seed(const string& s, int argc, char **argv) {
	const string	str_value = flag_value(s, argc, argv);
	if(str_value.empty()) {
		return -1;
	}
	else {
		const int	seed = std::stoi(str_value);
		if(seed > 0)
			return seed;
		else
			throw std::invalid_argument("seed must be positive");
	}
}

void Option::usage() {
	cerr << "usage : BreedingSim num_inds num_chroms [-s seed(> 0)]"
			<< " [-t num threads] [-o path out] [-p path pheno out]." << endl;
}

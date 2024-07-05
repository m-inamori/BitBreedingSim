#include <iostream>
#include <stdexcept>
#include "../include/option.h"

using namespace std;


//////////////////// Option ////////////////////

const Option *Option::create(int argc, char **argv) {
	if(argc != 3)
		return nullptr;
	
	try {
		const size_t	ni = std::stoull(argv[1]);
		const size_t	nc = std::stoull(argv[2]);
		return new Option(ni, nc);
	}
	catch(const invalid_argument& e) {
		return nullptr;
	}
}

void Option::usage() {
	cerr << "usage : BreedingSim num_inds num_chroms." << endl;
}

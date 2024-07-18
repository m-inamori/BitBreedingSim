#include <iostream>
#include <fstream>
#include <chrono>
#include "../include/population.h"
#include "../include/Map.h"
#include "../include/option.h"
#include "../include/common.h"

using namespace std;

void BreedingSim(const Option *op) {
	const Map *gmap = Map::create_default(op->num_chroms, 1000, 1.0);
	const auto	mat_origins = Population::create_origins(10, *gmap, "mat_");
	const auto	pat_origins = Population::create_origins(10, *gmap, "pat_");
auto start = std::chrono::high_resolution_clock::now();
	const auto	progs = Population::cross(op->num_inds, *mat_origins,
											*pat_origins, *gmap, "progs_",
											op->num_threads);
auto end = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> diff = end-start;
std::cout << "Time to run function: " << diff.count() << " s\n";
	if(!op->path_out.empty()) {
		ofstream	ofs(op->path_out.c_str());
		if(ofs) {
			progs->write(ofs);
		}
	}
	
	delete mat_origins;
	delete pat_origins;
	delete progs;
	delete gmap;
}

int main(int argc, char **argv) {
	const Option	*option = Option::create(argc, argv);
	if(option == nullptr) {
		Option::usage();
		exit(1);
	}
	
	BreedingSim(option);
	delete option;
	return 0;
}

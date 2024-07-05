#include <iostream>
#include <chrono>
#include "../include/individual.h"
#include "../include/option.h"

using namespace std;

void BreedingSim(const Option *op) {
	const auto	mat_origins = Individual::create_origins(10, op->num_chroms,
															1.0, 1000, "mat_");
	const auto	pat_origins = Individual::create_origins(10, op->num_chroms,
															1.0, 1000, "pat_");
auto start = std::chrono::high_resolution_clock::now();
	const auto	progs = Individual::cross(mat_origins,
											pat_origins, op->num_inds,
											"progs");
auto end = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> diff = end-start;
std::cout << "Time to run function: " << diff.count() << " s\n";
    cout << progs.size() << endl;
}

int main(int argc, char **argv) {
	const Option	*option = Option::create(argc, argv);
	if(option == nullptr) {
		Option::usage();
		exit(1);
	}
	
	BreedingSim(option);
	return 0;
}

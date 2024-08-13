#include <iostream>
#include <fstream>
#include <chrono>
#include "../include/BaseInfo.h"
#include "../include/population.h"
#include "../include/Map.h"
#include "../include/option.h"
#include "../include/common.h"

using namespace std;

void BreedingSim(const Option *op) {
	BaseInfo	*info = BaseInfo::create_default(op->seed);
	info->set_trait_AD_multi(10, 0.6, 0.7);
	
	const auto	mat_origins = Population::create_origins(10, info, "mat_");
	const auto	pat_origins = Population::create_origins(10, info, "pat_");
auto start = std::chrono::high_resolution_clock::now();
	const auto	progs = Population::cross(op->num_inds, *mat_origins,
											*pat_origins, info, "progs_",
											op->num_threads);
auto end = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> diff = end-start;
std::cout << "Time to run function: " << diff.count() << " s\n";
	
	// 表現型で選抜する
	const size_t	num = 10;
	vector<pair<double, size_t>>	sorted_phenos(progs->num_inds());
	const auto	phenos = progs->get_phenotypes(0);
	for(size_t i = 0; i < progs->num_inds(); ++i) {
		sorted_phenos[i] = make_pair(phenos[i], i);
	}
	std::sort(sorted_phenos.begin(), sorted_phenos.end());
	vector<size_t>	selected_indices(num);
	for(size_t i = 0; i < num; ++i) {
		selected_indices[i] = sorted_phenos[i+progs->num_inds()-num].second;
	}
	const auto	selected_progs = progs->select(selected_indices);
	cout << "selected mean : " << selected_progs->mean(0)
				<< " sd : " << selected_progs->stddev(0) << endl;
	
	// もう一度交配する
	const auto	progs2 = Population::cross(op->num_inds, *selected_progs,
											*selected_progs, info, "progs2_",
											op->num_threads);
	
#if 0
	if(!op->path_out.empty()) {
		ofstream	ofs(op->path_out.c_str());
		if(ofs) {
			progs2->write(ofs);
		}
	}
#else
	if(!op->path_out.empty()) {
		ofstream	ofs(op->path_out.c_str());
		if(ofs) {
			progs2->write_phenotypes(ofs);
		}
	}
#endif
	cout << "1st mean : " << progs->mean(0)
				<< " sd : " << progs->stddev(0) << endl;
	cout << "2nd mean : " << progs2->mean(0)
				<< " sd : " << progs2->stddev(0) << endl;
	
	mat_origins->dispay_QTLs(0);
	pat_origins->dispay_QTLs(0);
	progs->dispay_QTLs(0);
	selected_progs->dispay_QTLs(0);
	progs2->dispay_QTLs(0);
	
	delete mat_origins;
	delete pat_origins;
	delete progs;
	delete selected_progs;
	delete progs2;
	delete info;
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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "../include/individual.h"
#include "../include/Map.h"
#include "../include/common.h"

using namespace std;

Individual::~Individual() {
	Common::delete_all(haps1);
	Common::delete_all(haps2);
}

vector<const Individual *> Individual::cross(
							const vector<const Individual *>& mothers,
							const vector<const Individual *>& fathers,
							size_t num_inds, const string& name_base, int T) {
	vector<const Individual *>	inds(num_inds);
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(i, T, mothers, fathers,
											num_inds, name_base, inds);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
			(void *(*)(void *))&cross_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		cross_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	
	return inds;
}

void Individual::cross_in_thread(void *config) {
	auto	*c = (ConfigThread *)config;
	
	std::random_device	seed_gen;
	std::default_random_engine	engine;
cerr << "a " << c->first << endl;
	
	std::uniform_int_distribution<size_t>	dist1(0, c->mothers.size()-1);
	std::uniform_int_distribution<size_t>	dist2(0, c->fathers.size()-1);
cerr << "b " << c->first << endl;
	
	for(size_t i = c->first; i < c->num_inds; i += c->num_threads) {
		stringstream	ss;
		ss << setfill('0') << setw(8) << c->name_base << i + 1;
		const size_t	mother_index = dist1(engine);
		const size_t	father_index = dist2(engine);
		c->results[i] = cross_each(ss.str(), c->mothers[mother_index],
										c->fathers[father_index], engine);
	}
}

const Individual *Individual::cross_each(const string& name,
										 const Individual *mother,
										 const Individual *father,
										 std::default_random_engine &engine) {
	vector<const Haplotype *>	haps1;
	vector<const Haplotype *>	haps2;
	for(size_t i = 0; i < mother->num_chroms(); ++i) {
		const Haplotype	*hap1 = Haplotype::reduce(mother->haps1[i],
													mother->haps2[i], engine);
		haps1.push_back(hap1);
		const Haplotype	*hap2 = Haplotype::reduce(father->haps1[i],
													father->haps2[i], engine);
		haps2.push_back(hap2);
	}
	return new Individual(name, mother, father, haps1, haps2, mother->map);
}

vector<const Individual *> Individual::create_origins(
									size_t num_inds, size_t num_chroms,
									double length, size_t num_markers,
									const string& name_base) {
	// とりあえずデフォルトの地図
	shared_ptr<const Map>	map_(Map::create_default(num_chroms,
														num_markers, length));
	vector<const Haplotype *>	all_haps(num_inds * 2 * num_chroms);
	for(size_t i = 0; i < num_chroms; ++i) {
		const auto	haps = Haplotype::create_origin_haplotypes(num_inds*2,
														length, num_markers,
														*(map_->get_chr(i)));
		for(size_t j = 0; j < num_inds*2; ++j)
			all_haps[j*num_chroms+i] = haps[j];
	}
	
	vector<const Individual *>	inds;
	for(size_t i = 0; i < num_inds; ++i) {
		stringstream	ss;
		ss << setfill('0') << setw(3) << name_base << i + 1;
		vector<const Haplotype *>	haps1(all_haps.begin()+ i*2*num_chroms,
										  all_haps.begin()+(i*2+1)*num_chroms);
		vector<const Haplotype *>	haps2(all_haps.begin()+(i*2+1)*num_chroms,
										  all_haps.begin()+(i*2+2)*num_chroms);
		inds.push_back(new Individual(ss.str(), nullptr, nullptr,
												haps1, haps2, map_));
	}
	return inds;
}

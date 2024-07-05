#ifndef __INDIVIDUAL
#define __INDIVIDUAL

#include <random>
#include "haplotype.h"

class Map;


//////////////////// Individual ////////////////////

class Individual {
	const std::string	name;
	const Individual	*mother;
	const Individual	*father;
	const std::vector<const Haplotype *>	haps1;
	const std::vector<const Haplotype *>	haps2;
	std::shared_ptr<const Map>	map;
	
public:
	Individual(const std::string& name_, const Individual *m,
				const Individual *f,
				const std::vector<const Haplotype *>& h1,
				const std::vector<const Haplotype *>& h2,
				std::shared_ptr<const Map> m_) :
			name(name_), mother(m), father(f), haps1(h1), haps2(h2), map(m_) { }
	
	std::size_t num_chroms() const { return haps1.size(); }
	
	const Haplotype *reduce(std::default_random_engine &engine) const;
	
	static std::vector<const Individual *> cross(
								const std::vector<const Individual *>& mothers,
								const std::vector<const Individual *>& fathers,
								std::size_t num_inds,
								const std::string& name_base);
	static const Individual *cross_each(const std::string& name,
										const Individual *mother,
										const Individual *father,
										std::default_random_engine &engine);
	static std::vector<const Individual *> create_origins(
								std::size_t num_inds, std::size_t num_chroms,
								double length, std::size_t num_markers,
								const std::string& name_base);
};

#endif

#ifndef __OPTION
#define __OPTION


//////////////////// Option ////////////////////

class Option {
public:
	const std::size_t	num_inds;
	const std::size_t	num_chroms;
	
public:
	Option(std::size_t ni, std::size_t nc) : num_inds(ni), num_chroms(nc) { }
	
	static const Option *create(int argc, char **argv);
	static void usage();
};

#endif

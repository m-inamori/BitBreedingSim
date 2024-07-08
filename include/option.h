#ifndef __OPTION
#define __OPTION


//////////////////// Option ////////////////////

class Option {
public:
	const std::size_t	num_inds;
	const std::size_t	num_chroms;
	const int			num_threads;
	
public:
	Option(std::size_t ni, std::size_t nc, int T) :
						num_inds(ni), num_chroms(nc), num_threads(T) { }
	
	static const Option *create(int argc, char **argv);
	static std::string flag_value(const std::string& s, int argc, char **argv);
	static void usage();
};

#endif

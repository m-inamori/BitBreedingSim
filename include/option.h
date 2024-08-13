#ifndef __OPTION
#define __OPTION


//////////////////// Option ////////////////////

class Option {
public:
	const std::size_t	num_inds;
	const std::size_t	num_chroms;
	const int			seed;
	const int			num_threads;
	const std::string	path_out;
	
public:
	Option(std::size_t ni, std::size_t nc, int s, int T, const std::string& o) :
										num_inds(ni), num_chroms(nc),
										seed(s), num_threads(T), path_out(o) { }
	
	static const Option *create(int argc, char **argv);
	static std::string flag_value(const std::string& s, int argc, char **argv);
	static int flag_value_int(const std::string& s, int argc, char **argv,
															int default_value);
	static int flag_value_seed(const std::string& s, int argc, char **argv);
	static void usage();
};

#endif

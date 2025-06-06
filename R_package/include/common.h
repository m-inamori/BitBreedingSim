// common.h
#ifndef __COMMON
#define __COMMON

#include <vector>
#include <string>
#include <algorithm>


//////////////////// Common ////////////////////

namespace Common {
	std::string strip(const std::string& s);
	bool empty_line(const std::string& s);
	std::vector<std::string> split(const std::string& s, char delim);
	std::string join(const std::vector<std::string>& v, char delim);
	void chop(char *buff);
	std::vector<std::string> merge_vector(const std::vector<std::string>& v1,
										  const std::vector<std::string>& v2);
	
	std::vector<std::vector<std::string>> read_csv(
							const std::string& path, char delim=',');
	std::vector<std::vector<std::string>> read_tsv(const std::string& path);
	// read whitespace separeted file
	std::vector<std::vector<std::string>> read_wsv(const std::string& path);
	void write_tsv(const std::vector<std::string>& v, std::ostream& os);
	
	// e.g. [2, 5, 7], [3, 5, 6] -> [2, 5, 7, 3, 6]
	std::vector<std::string> merge_vectors(
								const std::vector<std::string>& v1,
								const std::vector<std::string>& v2);
	bool is_all_same(const std::string& seq);
	template<typename T>
	bool is_all_same(const std::vector<T>& v) {
		for(auto p = v.begin() + 1; p != v.end(); ++p) {
			if(*p != v.front())
				return false;
		}
		return true;
	}
	
	template<typename T>
	std::vector<T> unique_vector(const std::vector<T>& v) {
		std::vector<T>	w(1U, v.front());
		for(auto p = v.begin() + 1; p != v.end(); ++p) {
			if(std::find(w.begin(), w.end(), *p) == w.end())
				w.push_back(*p);
		}
		return w;
	}
	
	template<typename T>
	T dot_product(const std::vector<T>& v, const std::vector<T>& w) {
		T	s = 0;
		for(size_t i = 0; i < v.size(); ++i) {
			s += v[i] * w[i];
		}
		return s;
	}
	
	template<typename T>
	std::vector<T> multiply_by_constant(const std::vector<T>& v, T a) {
		std::vector<T>	w(v.size());
		for(size_t i = 0; i < v.size(); ++i) {
			w[i] = v[i] * a;
		}
		return w;
	}
	
	template<typename T>
	void connect_vector(const std::vector<T>& v, const std::vector<T>& w,
													std::vector<T>& dest) {
		dest.reserve(v.size() + w.size());
		dest.insert(dest.end(), v.begin(), v.end());
		dest.insert(dest.end(), w.begin(), w.end());
	}
	
	template<typename T>
	void delete_all(const std::vector<T *>& v) {
		for(auto p = v.begin(); p != v.end(); ++p)
			delete *p;
	}
	
	bool is_int(const std::string& s);
	bool is_double(const std::string& s);
}
#endif

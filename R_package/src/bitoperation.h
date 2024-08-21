#pragma once

#include "Int.h"

namespace BitOperation {
	// ex) 3 -> 1...111000
	inline Int::ull upper_mask(std::size_t k) {
		if(k == 0)
			return ~0;
		else
			return ((1ULL << (64 - k)) - 1) << k;
	}
	
	// ex) 3 -> 0...000111
	inline Int::ull lower_mask(std::size_t k) {
		return (1ULL << k) - 1;
	}
	
	// ex) 3, 5 -> 0...00011000
	inline Int::ull middle_mask(std::size_t first, std::size_t last) {
		return ((1ULL << (last - first)) - 1) << first;
	}
}

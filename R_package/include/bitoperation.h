#pragma once

namespace BitOperation {
	// ex) 3 -> 1...111000
	inline uint64_t upper_mask(std::size_t k) {
		if(k == 0)
			return ~0;
		else
			return ((1ULL << (64 - k)) - 1) << k;
	}
	
	// ex) 3 -> 0...000111
	inline uint64_t lower_mask(std::size_t k) {
		return (1ULL << k) - 1;
	}
	
	// ex) 3, 5 -> 0...00011000
	inline uint64_t middle_mask(std::size_t first, std::size_t last) {
		return ((1ULL << (last - first)) - 1) << first;
	}
}

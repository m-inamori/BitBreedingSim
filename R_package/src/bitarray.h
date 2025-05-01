#pragma once

#include <vector>
#include <cstdint>

namespace BitArray {
	using ConstIter = std::vector<uint64_t>::const_iterator;
	using Iter = std::vector<uint64_t>::iterator;
	
	inline int get(ConstIter source, std::size_t i) {
		const std::size_t	q = i >> 6;
		const std::size_t	r = i & 63;
		return static_cast<int>((*(source + q) >> r) & 1);
	}
	void copy(ConstIter source, std::size_t first, std::size_t last, Iter dest);
}

#include "../include/bitarray.h"
#include "../include/bitoperation.h"

using namespace std;

void BitArray::copy(ConstIter source, size_t first, size_t last, Iter dest) {
	const size_t	first_q = first / 64;
	const size_t	first_r = first % 64;
	const size_t	last_q = (last - 1) / 64;
	
	const Int::ull	mask = BitOperation::upper_mask(first_r);
	*(dest + first_q) &= (~0) ^ mask;
	*(dest + first_q) |= (*(source + first_q)) & mask;
	if(first_q < last_q) {
		std::copy(source + first_q + 1, source + last_q + 1,
											dest + first_q + 1);
	}
}

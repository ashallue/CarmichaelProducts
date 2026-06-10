/* I got this idea from Jon Grantham, though I don't know if it 
orginates with him.  The idea is rather than storing a prime, to instead store the 
bit representations of the exponents of p-1 as a divisor of Lambda.  The advantage
is a significant savings in space.  For the example of the billion prime case, primes 
go from 160 bit ZZ's (probably take even more than that) to 64 bit longs.  The 
disadvantages are first, a significant slowdown in converting back and forth and 
second, the inability to represent an arbitrary integer modulo Lambda in this fashion.
In particular, this means we can't multiply CondensedInteger's.

This code was written by Steven Hayman and Andrew Shallue in 2011, then in Summer 2026 modified and updated 
from ZZ to mpz_t by Andrew Shallue with the help of Gemini Flash-3.5
*/

#include <vector>
#include <gmp.h>
#include <math.h>

#ifndef CONDENSEDINTEGER_H
#define CONDENSEDINTEGER_H

class CondensedInteger {

	private:
		std::vector<bool> rep;

		// lambda represented as the product of primes to exponents
		static std::vector<long> primes;
		static std::vector<long> exponents;

	public:
		// decided against pass-by-ref here. These are comparatively small, and copies are safer
		static void set(std::vector<long> primes_, std::vector<long> exponents_);

	CondensedInteger();
  
	// Accept the mpz_t by const reference/pointer (Gemini claims GMP handles this internally)
	// main action here is computing p-1 as a divisor of lambda, then storing bits of the exponents in rep
	CondensedInteger(const mpz_t p);

	// rop is the destination variable.  Again, apparently pass-by-ref handled internally
	void to_mpz(mpz_t rop);
};

#endif